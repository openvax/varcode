"""Bounded, sequence-only discovery of starts created in paired RNA products.

This consumes existing phase/splice products. It does not infer initiation,
expression, downstream regulation, or protein absence in the baseline.
"""

from collections import defaultdict
from dataclasses import dataclass
from hashlib import sha256
import json

from serializable import DataclassSerializable

from .effect_hypotheses import RealizedTranscriptProduct
from .effect_traversal import candidate_children, iter_effects
from .effects.codon_tables import codon_table_for_transcript
from .effects.effect_classes import MutationEffect, StartGain
from .mutant_transcript import MutantTranscript
from .orf import _plain
from .reference_metadata import reference_metadata


@dataclass(frozen=True)
class ORFDiscoveryPolicy(DataclassSerializable):
    """Explicit limits and candidate initiators; no Kozak exclusion threshold.

    ``start_codons`` are hypotheses, including caller-selected non-AUG codons.
    Their first residue is conditionally methionine. Default discovery is AUG
    only, independently of a genetic code's larger set of known initiators.
    Limits apply per paired product, except ``max_products`` per input effect.
    """

    start_codons: tuple = ("ATG",)
    max_scan_bases: int = 100_000
    max_protein_length: int = 1_000
    max_candidates: int = 32
    max_products: int = 256

    def __post_init__(self):
        codons = tuple(self.start_codons)
        if not codons or len(set(codons)) != len(codons) or any(
                not isinstance(c, str) or len(c) != 3 or set(c) - set("ACGT") for c in codons):
            raise ValueError("start_codons must be distinct, unambiguous DNA triplets")
        object.__setattr__(self, "start_codons", codons)
        for name in ("max_scan_bases", "max_protein_length", "max_candidates", "max_products"):
            value = getattr(self, name)
            if type(value) is not int or value < 1:
                raise ValueError("%s must be a positive integer" % name)


@dataclass(frozen=True)
class DiscoveredORF(DataclassSerializable):
    """A mutant RNA coding hypothesis, not a reference ORF or an observation.

    RNA intervals are zero-based, half-open; ``end`` includes a terminal stop
    when present. Base origins use the caller's reference coordinates.
    """

    orf_id: str
    transcript_id: str
    reference: dict
    start: int
    end: int
    start_codon: str
    translation_table: int
    termination: str
    stop_offset: int | None
    base_origins: tuple
    initiation_context: dict
    baseline_codons: tuple
    creation: str
    source: str = "varcode_created_start"
    orf_type: str = "variant_created"
    kind: str = "discovered"


@dataclass(frozen=True)
class ORFDiscoveryAssessment(DataclassSerializable):
    """One paired comparison, including empty, bounded and unavailable searches."""

    status: str
    reasons: tuple
    effects: tuple
    scanned_bases: int
    baseline_cdna_sequence: str | None
    mutant_cdna_sequence: str | None
    source_evidence: dict


@dataclass(frozen=True)
class ORFDiscovery(DataclassSerializable):
    """Discovery results and their declared scope; an empty set is not no ORFs."""

    policy: ORFDiscoveryPolicy
    assessments: tuple

    @property
    def effects(self):
        return tuple(effect for item in self.assessments for effect in item.effects)

    @property
    def complete(self):
        """Whether every supplied product was fully searched under this policy."""
        return bool(self.assessments) and all(item.status == "complete" for item in self.assessments)


def _origins(product):
    """Normalize explicit ancestry, never infer it by sequence alignment."""
    evidence = product.evidence
    origins = evidence.get("base_origins")
    if origins is not None:
        result = []
        for origin in origins:
            if origin is None:
                result.append(None)
            elif (len(origin) == 3 and isinstance(origin[0], str)
                  and type(origin[1]) is int and origin[1] > 0 and origin[2] in ("+", "-")):
                result.append(tuple(origin))
            else:
                raise ValueError("base_origins require (contig, positive position, strand) or None")
        return tuple(result)
    positions = evidence.get("genomic_positions")
    transcript = product.transcript
    if positions is None or transcript is None:
        return None
    strand = "-" if transcript.on_backward_strand else "+"
    return tuple(None if p is None else (str(transcript.contig), p, strand) for p in positions)


def _translate(sequence, start, table, limit):
    amino_acids = []
    for offset in range(start, len(sequence) - 2, 3):
        codon = sequence[offset:offset + 3]
        if codon in table.stop_codons:
            return "".join(amino_acids), offset + 3, offset, "observed_stop"
        if len(amino_acids) == limit:
            return "".join(amino_acids), offset, None, "protein_length_limit"
        if codon not in table.forward_table:
            return "".join(amino_acids), offset, None, "ambiguous_sequence"
        amino_acids.append("M" if offset == start else table.forward_table[codon])
    return "".join(amino_acids), start + 3 * len(amino_acids), None, "sequence_boundary"


def _context(sequence, start):
    return dict(upstream=sequence[max(0, start - 6):start],
                start_codon=sequence[start:start + 3], downstream=sequence[start + 3:start + 9],
                upstream_bases_available=min(start, 6),
                biological_five_prime_completeness="unknown",
                interpretation="sequence_prior_only; initiation_not_measured")


def discover_created_orfs_in_products(baseline, mutant, *, variant=None, policy=None,
                                      hypothesis=None, source_evidence=None):
    """Search a paired, oriented RNA comparison and return one assessment.

    Both inputs are ``RealizedTranscriptProduct`` objects for the same host
    transcript/reference. ``evidence.base_origins`` contains one mapped
    ``(contig, position, strand)`` or None per base. Shared transcript-model
    ``genomic_positions`` are also supported. None means unknown ancestry unless
    its offset is explicitly in ``inserted_base_offsets``. Inserted starts with
    inserted sequence in the baseline remain ambiguous without shared edit IDs.

    Partial baseline RNA cannot establish a newly included start; it is reported
    unavailable. Complete supplied RNA does not establish the biological 5' end.
    Sequence-only candidates retain the source hypothesis and paired RNA.
    """
    policy = ORFDiscoveryPolicy() if policy is None else policy
    if not isinstance(policy, ORFDiscoveryPolicy):
        raise TypeError("policy must be ORFDiscoveryPolicy")
    if not isinstance(baseline, RealizedTranscriptProduct) or not isinstance(mutant, RealizedTranscriptProduct):
        raise TypeError("paired inputs must be RealizedTranscriptProduct objects")
    before, after = baseline.cdna_sequence, mutant.cdna_sequence
    evidence = dict(source_evidence or {})
    evidence.update(hypothesis=_plain(hypothesis), policy=policy.to_dict(),
                    baseline_evidence=_plain(baseline.evidence), mutant_evidence=_plain(mutant.evidence),
                    mutant_translation_evidence="not_assessed", downstream_cds_regulation="not_assessed")

    def assessment(status, reasons, effects=(), scanned=0):
        return ORFDiscoveryAssessment(status, tuple(reasons), tuple(effects), scanned, before, after, evidence)

    if before is None or after is None:
        return assessment("unavailable", ("paired_rna_sequence_unavailable",))
    transcript = mutant.transcript
    if transcript is None or baseline.transcript is None:
        return assessment("unavailable", ("host_transcript_unavailable",))
    reference = reference_metadata(transcript.genome)
    if (baseline.transcript.id != transcript.id
            or reference_metadata(baseline.transcript.genome) != reference):
        raise ValueError("Paired products must use the same host transcript and reference identity")
    if any(product.evidence.get("partial") or product.evidence.get("partial_sequence")
           for product in (baseline, mutant)):
        return assessment("unavailable", ("partial_rna_comparison",))
    if len(before) > policy.max_scan_bases:
        return assessment("bounded", ("baseline_scan_limit",))
    old_origins, new_origins = _origins(baseline), _origins(mutant)
    if old_origins is None or new_origins is None:
        return assessment("unavailable", ("base_origins_unavailable",))
    if len(old_origins) != len(before) or len(new_origins) != len(after):
        raise ValueError("A base origin is required for every RNA nucleotide")
    before, after = before.upper(), after.upper()
    table = codon_table_for_transcript(transcript)
    if set(policy.start_codons) & table.stop_codons:
        raise ValueError("A discovery start codon is a stop in the selected genetic code")
    triples = defaultdict(list)
    for i in range(len(before) - 2):
        key = old_origins[i:i + 3]
        if None not in key:
            triples[key].append(before[i:i + 3])
    inserted = set(mutant.evidence.get("inserted_base_offsets", ()))
    if any(type(i) is not int or not 0 <= i < len(after) or new_origins[i] is not None for i in inserted):
        raise ValueError("inserted_base_offsets must identify unmapped bases in mutant RNA")
    effects, reasons = [], []
    scanned = min(len(after), policy.max_scan_bases)
    if scanned < len(after):
        reasons.append("mutant_scan_limit")
    if set(before) - set("ACGT") or set(after[:scanned]) - set("ACGT"):
        reasons.append("ambiguous_rna_sequence")
    for start in range(max(0, scanned - 2)):
        codon = after[start:start + 3]
        if codon not in policy.start_codons:
            continue
        origins = new_origins[start:start + 3]
        unmapped = [start + j for j, origin in enumerate(origins) if origin is None]
        if unmapped and (not set(unmapped) <= inserted or None in old_origins):
            reasons.append("ambiguous_start_ancestry")
            continue
        old_codons = tuple(triples.get(origins, ())) if not unmapped else ()
        if any(c in policy.start_codons for c in old_codons):
            continue
        if any(set(c) - set("ACGT") for c in old_codons):
            reasons.append("ambiguous_baseline_codon")
            continue
        # With an incomplete baseline origin map, absence of the same triplet
        # cannot prove a new junction/included segment.
        if not old_codons and None in old_origins:
            reasons.append("ambiguous_baseline_ancestry")
            continue
        if len(effects) == policy.max_candidates:
            reasons.append("candidate_limit")
            break
        protein, end, stop, termination = _translate(after, start, table, policy.max_protein_length)
        creation = "changed_codon" if old_codons else "new_rna_adjacency_or_sequence"
        identity = dict(reference=reference, transcript_id=transcript.id,
                        mutant_rna=after, start=start, end=end,
                        origins=new_origins[start:end], translation_table=table.ncbi_id)
        identifier = "created-orf:" + sha256(json.dumps(identity, sort_keys=True).encode()).hexdigest()[:24]
        model = DiscoveredORF(identifier, transcript.id, reference, start, end, codon,
                              table.ncbi_id, termination, stop, new_origins[start:end],
                              _context(after, start), old_codons, creation)
        candidate_evidence = {**evidence, "discovered_orf": model.to_dict(),
                              "baseline_protein_status": "initiation_not_assessed"}
        baseline_model = MutantTranscript.from_sequence(
            before, reference_transcript=transcript, annotator_name="created_orf",
            evidence={**candidate_evidence, "base_origins": old_origins})
        mutant_model = MutantTranscript.from_sequence(
            after, reference_transcript=transcript, mutant_protein_sequence=protein,
            annotator_name="created_orf", evidence={**candidate_evidence, "base_origins": new_origins})
        effects.append(StartGain(variant, transcript, model, baseline_model, mutant_model,
                                 policy, hypothesis=hypothesis, source_evidence=evidence))
    reasons = tuple(dict.fromkeys(reasons))
    status = "bounded" if any(reason.endswith("_limit") for reason in reasons) else "partial" if reasons else "complete"
    return assessment(status, reasons, effects, scanned)


def discover_created_orfs(effect, *, policy=None):
    """Discover new starts in every existing paired outcome of one effect.

    The ordinary annotation is unchanged. Returned ``effects`` expose the common
    sequence/host fields plus ``discovered_orf``. Realized phase alternatives
    remain separate. Legacy/unmapped results report unavailable discovery.
    """
    if not isinstance(effect, MutationEffect):
        raise TypeError("discover_created_orfs expects one MutationEffect")
    policy = ORFDiscoveryPolicy() if policy is None else policy
    if not isinstance(policy, ORFDiscoveryPolicy):
        raise TypeError("policy must be ORFDiscoveryPolicy")
    assessments, seen = [], set()
    for node in iter_effects(effect):
        candidates = getattr(node, "realized_candidates", ())
        for candidate in candidates:
            for outcome in candidate.outcomes:
                if id(outcome) in seen:
                    continue
                seen.add(id(outcome))
                if len(assessments) == policy.max_products:
                    assessments.append(ORFDiscoveryAssessment(
                        "bounded", ("product_limit",), (), 0, None, None, {}))
                    return ORFDiscovery(policy, tuple(assessments))
                assessments.append(discover_created_orfs_in_products(
                    outcome.baseline, outcome.mutant, variant=node.variant, policy=policy,
                    hypothesis=outcome.hypothesis,
                    source_evidence=dict(source_effect=type(outcome.effect).__name__,
                                         probability=candidate.probability, phase_mass=candidate.phase_mass)))
        if not candidates and not candidate_children(node):
            assessments.append(ORFDiscoveryAssessment(
                "unavailable", ("paired_products_unavailable",), (), 0, None, None,
                dict(source_effect=type(node).__name__)))
    return ORFDiscovery(policy, tuple(assessments))
