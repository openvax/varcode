"""RNA-backed effect comparisons using the shared openvax-v2 osteosarc bundle.

Varcode predicts coding consequences from DNA variants, germline context, and
externally supplied RNA phasing/evidence (see ``docs/library_roles.md``);
BAM-level read phasing and RNA assembly live in Isovar. This suite verifies:

1. The neutral shared ``openvax-v2`` osteosarc bundle keeps complete templates
   and RNA source provenance, reuses the shared cache offline, and rejects
   missing or tampered bundle manifests on an isolated copy.
2. Independent NCBI RefSeq references (``tests/data/sid_rna_effect_references.json``)
   match the Ensembl 81 NTF3 and MAP2 transcripts.
3. Passing audited RNA phasing observations through Varcode's public
   ``ReadPhasingSource`` / ``MolecularPhaseResolver`` interface produces the
   independently audited compound proteins and cDNA edits, survives JSON
   round-trips, preserves every unphased DNA prediction (never treating
   unresolved or absent RNA as absence of a DNA event), and does not conflate
   competing MAP2 alleles.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection,
    MolecularPhaseResolver,
    ReadPhasingSource,
    Variant,
    VariantCollection,
)
from varcode.effects.effect_classes import HaplotypeEffect, Silent, Substitution


SHORT_SOURCE = "2024.06.11.bostongene.align.tcga.protocol.dr32Aligned.sorted"
LONG_SOURCE = "IPISRC044_T1_sclrs_ONT.tagged"
SHORT_NTF3 = SHORT_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_NTF3 = LONG_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_EXOC4 = LONG_SOURCE + ".EXOC4-chr7-133274996"
MAP2_PREFIX = (
    "isovar/osteosarc/figure_comparisons/corpus/"
    "context-evidence.json.gz#/map2/sources/"
)
SHORT_MAP2_T1 = MAP2_PREFIX + "11/original_records"
SHORT_MAP2_T2 = MAP2_PREFIX + "12/original_records"
LONG_MAP2_T2 = MAP2_PREFIX + "3/original_records"
MAP2_MEMBERS = (SHORT_MAP2_T1, SHORT_MAP2_T2, LONG_MAP2_T2)
MEMBERS = (SHORT_NTF3, LONG_NTF3, LONG_EXOC4) + MAP2_MEMBERS

MANIFEST_SHA256 = (
    "b42dbce529cb5724153f035fca3783a41e673b31610a30851b2ef96d4ba43734"
)
REFERENCE_PATH = Path(__file__).parent / "data" / "sid_rna_effect_references.json"
REFERENCE_SHA256 = (
    "d9e21cc96450d4a14e7b41d0900df48792c9f64b103b9ce01ab5cc841508d29e"
)

NTF3 = ("12", 5494381, "A", "G")
NTF3_COMPOUND_PARTNER = ("12", 5494382, "G", "T")
# Read-observed companion alleles; no germline/somatic status is assumed.
NTF3_COMPANION = ("12", 5494466, "G", "A")
EXOC4 = ("7", 133274996, "G", "T")
EXOC4_COMPANION = ("7", 133895694, "G", "A")
MAP2_PARTS = (
    ("2", 209694769, "C", "A"),
    ("2", 209694770, "T", "G"),
    ("2", 209694773, "GCTACTGTGTGTTCAATAAGTACACAGT", ""),
)
MAP2_COMPOUND = ("2", 209694768, "CCTGGGCTACTGTGTGTTCAATAAGTACACAGT", "CAGGG")
MAP2_HISTORICAL = ("2", 209694768, "CCTGGGCTACTGTGTGTTCAATA", "C")


class _AuditedRNAPhaseSource:
    """Minimal ``ReadPhasingSource`` delivering audited RNA phasing calls.

    Mirrors the contract Isovar's ``IsovarReadPhasing`` provides to
    ``MolecularPhaseResolver`` without implementing BAM parsing inside Varcode.
    """

    source = "rna_reads"

    def __init__(self, *, supported=(), cis_groups=(), trans_pairs=()):
        self._supported = set(supported)
        self._cis_groups = tuple(tuple(group) for group in cis_groups)
        self._trans_pairs = {frozenset(pair) for pair in trans_pairs}

    def has_evidence(self, variant) -> bool:
        return variant in self._supported

    def partners_in_cis(self, variant):
        if not self.has_evidence(variant):
            return ()
        partners = []
        for group in self._cis_groups:
            if variant in group:
                for member in group:
                    if member != variant and member not in partners:
                        partners.append(member)
        return tuple(partners)

    def in_cis(self, v1, v2, transcript=None):
        if v1 == v2:
            return True if self.has_evidence(v1) else None
        for group in self._cis_groups:
            if v1 in group and v2 in group:
                return True
        if frozenset((v1, v2)) in self._trans_pairs:
            return False
        return None


@pytest.fixture(scope="module")
def shared_rna_bundle(tmp_path_factory):
    import osteosarc
    import pysam
    from osteosarc.records import read_records, read_template

    root = Path(osteosarc.fetch_bundle("openvax-v2"))
    manifest_bytes = (root / "manifest.json").read_bytes()
    assert hashlib.sha256(manifest_bytes).hexdigest() == MANIFEST_SHA256
    recipe = json.loads((root / "recipe.json").read_text())
    assert recipe["selection"]["supplementary"] is True
    compound = recipe["targets"]["NTF3-chr12-5494381-compound"]
    assert (
        compound["assembly"],
        compound["contig"],
        compound["position"],
        compound["ref"],
        compound["alt"],
    ) == ("GRCh38", "chr12", 5494381, "AG", "GT")

    destination = tmp_path_factory.mktemp("sid-rna-effects")
    paths = osteosarc.export_bundle("openvax-v2", destination, members=MEMBERS)
    manifest = json.loads(manifest_bytes)
    # The legacy MAP2 members are regional acquisitions. Include every other
    # original record of their selected templates available in the bundle,
    # including distant mates and supplementary placements.
    for index, member in enumerate(MAP2_MEMBERS):
        selected = {record.template for record in read_records(paths[member])}
        source = manifest["members"][member]["source"]
        complete_path = destination / ("complete-map2-%d.bam" % index)
        with pysam.AlignmentFile(root / manifest["sources"][source]["bam"], "rb") as bam:
            with pysam.AlignmentFile(complete_path, "wb", template=bam) as output:
                for read in bam:
                    if read_template(read) in selected:
                        output.write(read)
        pysam.index(str(complete_path))
        paths[member] = complete_path
    return paths, root, manifest, recipe


@pytest.fixture(scope="module")
def effect_references():
    """Offline RefSeq references and separately audited protein expectations."""
    payload = REFERENCE_PATH.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == REFERENCE_SHA256
    return json.loads(payload)


def transcript_effects(effects, transcript_id, *, joint=False):
    return [
        effect
        for effect in effects
        if getattr(getattr(effect, "transcript", None), "id", None) == transcript_id
        and (len(getattr(effect, "variants", ())) > 1) is joint
    ]


def individual_signatures(effects):
    """Snapshot all individual predictions, including other transcript models."""
    return sorted(
        (
            str(effect.variant),
            getattr(getattr(effect, "transcript", None), "id", ""),
            type(effect).__name__,
            effect.short_description,
            effect.mutant_protein_sequence,
        )
        for effect in effects
        if len(getattr(effect, "variants", ())) < 2
    )


def assert_primary_protein(effect, case, references):
    expected = references["expectations"][case]
    reference = references["references"][expected["gene"]]
    protein = (
        reference["protein"][:expected["protein_prefix_length"]]
        + expected["protein_replacement"]
        + reference["protein"][expected["protein_suffix_start"]:]
    )
    assert len(protein) == expected["protein_length"]
    assert hashlib.sha256(protein.encode()).hexdigest() == expected["protein_sha256"]
    if case == "NTF3_companion":
        # Silent effects intentionally expose no mutant protein sequence.
        assert isinstance(effect, Silent)
        assert effect.aa_pos == 83 and effect.aa_ref == "P"
        assert effect.modifies_protein_sequence is False
        return
    assert effect.mutant_protein_sequence == protein


def assert_primary_cdna(effect, case, references):
    """Verify an available joint sequence product against the RefSeq edit audit."""
    expected = references["expectations"][case]
    reference = references["references"][expected["gene"]]
    if effect.mutant_transcript is not None:
        cdnas = [effect.mutant_transcript.cdna_sequence]
    else:
        cdnas = [
            outcome.mutant.cdna_sequence
            for candidate in effect.realized_candidates
            for outcome in candidate.outcomes
        ]
    assert cdnas
    start = min(effect.transcript.start_codon_spliced_offsets)
    delta = sum(len(alt) - len(ref) for _, ref, alt in expected["cds_edits"])
    for cdna in cdnas:
        cds = cdna[start:start + len(reference["cds"]) + delta]
        assert hashlib.sha256(cds.encode()).hexdigest() == expected["mutant_cds_sha256"]
        assert cds[expected["aa_pos"] * 3:expected["aa_pos"] * 3 + 3] == expected["alt_codon"]
        stop = expected["stop_codon_offset"]
        assert cds[stop:stop + 3] == expected["stop_codon"]


@pytest.mark.parametrize("member,records,templates,supplementary", [
    (SHORT_NTF3, 12, 6, 0),
    (LONG_NTF3, 15, 15, 0),
    (LONG_EXOC4, 40, 38, 2),
    (SHORT_MAP2_T1, 125, 78, 0),
    (SHORT_MAP2_T2, 26, 18, 0),
    (LONG_MAP2_T2, 33, 32, 1),
])
def test_shared_rna_bundle_members_keep_complete_templates_and_rna_provenance(
    shared_rna_bundle, member, records, templates, supplementary
):
    from osteosarc.records import read_records

    paths, root, manifest, _ = shared_rna_bundle
    source = manifest["members"][member]["source"]
    source_meta = manifest["sources"][source]
    assert source_meta["identity"]["assembly"] == "GRCh38"
    source_url = source_meta["identity"]["identity"]["url"]
    assert source_url.startswith(
        "https://sid-sijbrandij-osteosarc-dataset.s3.us-west-2.amazonaws.com/"
    )
    assert any(token in source_url for token in ("/RNA/", "/rna/", "/ONT/"))

    exported = list(read_records(paths[member]))
    identities = Counter(record.digest for record in exported)
    member_identities = Counter(manifest["members"][member]["records"])
    assert not member_identities - identities
    if member not in MAP2_MEMBERS:
        assert identities == member_identities
    selected = {record.template for record in exported}
    source_records = list(read_records(root / source_meta["bam"]))
    assert selected == {
        record.template
        for record in source_records
        if record.digest in member_identities
    }
    complete = Counter(
        record.digest for record in source_records if record.template in selected
    )
    assert identities == complete
    assert len(exported) == records
    assert len(selected) == templates
    assert sum(record.read.is_supplementary for record in exported) == supplementary
    if member == SHORT_NTF3:
        for template in selected:
            assert {
                record.read.flag & 0xC0
                for record in exported
                if record.template == template
            } == {0x40, 0x80}


def test_shared_rna_bundle_offline_reuse_and_tamper_rejection(
    shared_rna_bundle, tmp_path
):
    import osteosarc
    from osteosarc import Cache, IntegrityError, OfflineError

    _, root, _, _ = shared_rna_bundle
    # Reopen the already-cached shared bundle strictly offline.
    reopened = Path(osteosarc.fetch_bundle("openvax-v2", offline=True))
    assert reopened == root

    # An empty cache must refuse network downloads in offline mode.
    empty_cache = Cache(tmp_path / "empty", offline=True)
    with pytest.raises(OfflineError):
        osteosarc.fetch_bundle("openvax-v2", cache=empty_cache, offline=True)

    # Tamper only with a private copy, never the user's shared cache.
    isolated_cache = Cache(tmp_path / "tampered", offline=True)
    isolated_bundle = isolated_cache.workspace / "bundles" / root.name
    shutil.copytree(root, isolated_bundle)
    manifest_path = isolated_bundle / "manifest.json"
    manifest_path.chmod(0o644)
    tampered = bytearray(manifest_path.read_bytes())
    tampered[0] ^= 1
    manifest_path.write_bytes(tampered)
    with pytest.raises(IntegrityError):
        osteosarc.fetch_bundle("openvax-v2", cache=isolated_cache, offline=True)


@pytest.mark.parametrize("gene", ["NTF3", "MAP2"])
def test_primary_effect_references_match_the_selected_transcripts(
    effect_references, gene
):
    reference = effect_references["references"][gene]
    transcript = cached_release(81).transcript_by_id(reference["transcript_id"])
    assert transcript.coding_sequence == reference["cds"]
    assert transcript.protein_sequence == reference["protein"]
    assert (
        hashlib.sha256(transcript.coding_sequence.encode()).hexdigest()
        == reference["cds_sha256"]
    )
    assert (
        hashlib.sha256(transcript.protein_sequence.encode()).hexdigest()
        == reference["protein_sha256"]
    )
    window = reference["genomic_window"]
    start = min(transcript.start_codon_spliced_offsets)
    assert transcript.spliced_offset(window["start"]) - start == window["cds_offset"]
    assert (
        transcript.sequence[
            transcript.spliced_offset(window["start"]):
            transcript.spliced_offset(window["end"]) + 1
        ]
        == window["sequence"]
    )


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_ntf3_cis_rna_phase_predicts_compound_codon(effect_references, annotator):
    variants = tuple(
        Variant(*allele, genome=cached_release(81))
        for allele in (NTF3, NTF3_COMPOUND_PARTNER)
    )
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = _AuditedRNAPhaseSource(supported=variants, cis_groups=(variants,))
    assert isinstance(source, ReadPhasingSource)
    resolver = MolecularPhaseResolver(source)
    assert resolver.in_cis(*variants) is True
    assert resolver.in_cis(*reversed(variants)) is True
    assert tuple(resolver.phased_partners(variants[0])) == (variants[1],)
    assert tuple(resolver.phased_partners(variants[1])) == (variants[0],)

    effects = collection.effects(annotator=annotator, phase_resolver=resolver)
    assert individual_signatures(effects) == baseline

    transcript_id = effect_references["references"]["NTF3"]["transcript_id"]
    (joint,) = transcript_effects(effects, transcript_id, joint=True)
    assert isinstance(
        joint, Substitution if annotator == "transcript_model" else HaplotypeEffect
    )
    assert joint.variants == variants
    assert joint.phase_source == "rna_reads"
    assert joint.annotator == effects.annotator == annotator
    assert joint.annotator_version == effects.annotator_version
    assert_primary_protein(joint, "NTF3_compound", effect_references)
    assert_primary_cdna(joint, "NTF3_compound", effect_references)

    individuals = {
        effect.variant: effect
        for effect in transcript_effects(effects, transcript_id)
    }
    for variant, case in zip(variants, ("NTF3_first", "NTF3_second")):
        assert_primary_protein(individuals[variant], case, effect_references)

    restored = EffectCollection.from_json(effects.to_json())
    (restored_joint,) = transcript_effects(restored, transcript_id, joint=True)
    assert restored_joint.variants == variants
    assert restored_joint.phase_source == joint.phase_source
    assert_primary_protein(restored_joint, "NTF3_compound", effect_references)
    if annotator == "transcript_model":
        from .test_transcript_model_serialization import (
            assert_realized_annotation_equal,
        )

        assert_realized_annotation_equal(joint, restored_joint)
        assert_primary_cdna(restored_joint, "NTF3_compound", effect_references)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize(
    "alleles,cases,phase",
    [
        # Audited ONT NTF3 companion (0 both-alt, 9 first-only, 3 second-only, 2 both-ref).
        ((NTF3, NTF3_COMPANION), ("NTF3_first", "NTF3_companion"), False),
        # Insufficient fragment support (e.g. min_alt_reads=10 with 9 reads): unknown.
        ((NTF3, NTF3_COMPOUND_PARTNER), ("NTF3_first", "NTF3_second"), None),
    ],
)
def test_ntf3_trans_and_unknown_rna_phase_keep_separate_predictions(
    effect_references, annotator, alleles, cases, phase
):
    variants = tuple(Variant(*allele, genome=cached_release(81)) for allele in alleles)
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = _AuditedRNAPhaseSource(
        supported=variants if phase is False else (),
        trans_pairs=(variants,) if phase is False else (),
    )
    resolver = MolecularPhaseResolver(source)
    assert resolver.in_cis(*variants) is phase  # Unknown must remain distinct from trans.
    assert resolver.in_cis(*reversed(variants)) is phase
    assert all(not resolver.phased_partners(variant) for variant in variants)

    effects = collection.effects(annotator=annotator, phase_resolver=resolver)
    assert individual_signatures(effects) == baseline
    assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
    transcript_id = effect_references["references"]["NTF3"]["transcript_id"]
    individuals = {
        effect.variant: effect
        for effect in transcript_effects(effects, transcript_id)
    }
    assert set(individuals) == set(variants)
    for variant, case in zip(variants, cases):
        assert_primary_protein(individuals[variant], case, effect_references)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
def test_exoc4_unresolved_rna_phase_preserves_dna_predictions(annotator):
    # Audited ONT EXOC4 pair has individual expression (19 and 3 alt templates)
    # but only 1 shared alt fragment, leaving phase unresolved (None).
    variants = tuple(
        Variant(*allele, genome=cached_release(81))
        for allele in (EXOC4, EXOC4_COMPANION)
    )
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = _AuditedRNAPhaseSource(supported=variants)
    resolver = MolecularPhaseResolver(source)
    assert all(resolver.has_evidence(variant) for variant in variants)
    assert resolver.in_cis(*variants) is None
    assert resolver.in_cis(*reversed(variants)) is None

    effects = collection.effects(annotator=annotator, phase_resolver=resolver)
    assert individual_signatures(effects) == baseline
    assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
@pytest.mark.parametrize("phase", [True, None])
def test_map2_cis_and_unknown_rna_phase_for_complex_allele(
    effect_references, annotator, phase
):
    # Bulk T1 (4 templates) and T2 (2 templates) support the 3-part MAP2 allele
    # in cis; single ONT T2 witness (1 template) leaves phase unknown.
    variants = tuple(
        Variant(*allele, genome=cached_release(81)) for allele in MAP2_PARTS
    )
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = _AuditedRNAPhaseSource(
        supported=variants if phase is True else (),
        cis_groups=(variants,) if phase is True else (),
    )
    resolver = MolecularPhaseResolver(source)
    for left, right in ((0, 1), (0, 2), (1, 2)):
        assert resolver.in_cis(variants[left], variants[right]) is phase

    effects = collection.effects(annotator=annotator, phase_resolver=resolver)
    assert individual_signatures(effects) == baseline
    transcript_id = effect_references["references"]["MAP2"]["transcript_id"]
    individuals = {
        effect.variant: effect
        for effect in transcript_effects(effects, transcript_id)
    }
    assert set(individuals) == set(variants)
    for variant, case in zip(variants, ("MAP2_first", "MAP2_second", "MAP2_deletion")):
        assert_primary_protein(individuals[variant], case, effect_references)

    if phase is None:
        assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
        assert not any(resolver.has_evidence(variant) for variant in variants)
    else:
        (joint,) = transcript_effects(effects, transcript_id, joint=True)
        assert joint.variants == variants
        assert joint.phase_source == "rna_reads"
        assert joint.annotator == effects.annotator == annotator
        assert_primary_protein(joint, "MAP2_compound", effect_references)
        assert_primary_cdna(joint, "MAP2_compound", effect_references)
        (restored_joint,) = transcript_effects(
            EffectCollection.from_json(effects.to_json()), transcript_id, joint=True
        )
        assert restored_joint.variants == variants
        assert restored_joint.phase_source == joint.phase_source
        assert_primary_protein(restored_joint, "MAP2_compound", effect_references)
        assert_primary_cdna(restored_joint, "MAP2_compound", effect_references)


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
def test_map2_rna_support_does_not_conflate_competing_dna_alleles(
    effect_references, annotator
):
    variants = tuple(
        Variant(*allele, genome=cached_release(81))
        for allele in (MAP2_COMPOUND, MAP2_PARTS[-1], MAP2_HISTORICAL)
    )
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    # Only the corrected compound allele has RNA witnesses; zero matches for
    # the other two hypotheses must not merge or drop their DNA predictions.
    source = _AuditedRNAPhaseSource(supported=(variants[0],))
    resolver = MolecularPhaseResolver(source)
    assert resolver.has_evidence(variants[0]) is True
    assert resolver.has_evidence(variants[1]) is False
    assert resolver.has_evidence(variants[2]) is False
    assert resolver.in_cis(variants[0], variants[1]) is None
    assert resolver.in_cis(variants[0], variants[2]) is None

    effects = collection.effects(annotator=annotator, phase_resolver=resolver)
    assert individual_signatures(effects) == baseline
    assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
    transcript_id = effect_references["references"]["MAP2"]["transcript_id"]
    individuals = {
        effect.variant: effect
        for effect in transcript_effects(effects, transcript_id)
    }
    assert set(individuals) == set(variants)
    for variant, case in zip(
        variants, ("MAP2_compound", "MAP2_deletion", "MAP2_historical")
    ):
        assert_primary_protein(individuals[variant], case, effect_references)


def test_no_runtime_dependency_on_osteosarc_or_vaxrank():
    requirements = (
        Path(__file__).resolve().parents[1] / "requirements.txt"
    ).read_text()
    runtime_pkgs = {
        line.split("#", 1)[0].strip().split(">")[0].split("<")[0].split("=")[0].split("[")[0].strip().lower()
        for line in requirements.splitlines()
        if line.strip() and not line.strip().startswith("#")
    }
    assert "osteosarc" not in runtime_pkgs
    assert "vaxrank" not in runtime_pkgs
