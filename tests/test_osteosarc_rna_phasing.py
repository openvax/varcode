"""RNA phasing from complete templates in the published openvax-v2 bundle.

These are observed allele combinations, not independently validated tumor
haplotypes. See tests/README.md for the sources, manual fragment audit and
filtering policy. No alignments or bases are synthesized for these tests.
"""

from collections import Counter
import hashlib
import json
from pathlib import Path
from unittest.mock import Mock

import pytest
from pyensembl import cached_release

from varcode import (
    EffectCollection, MolecularPhaseResolver, RNAReadPhasingSource, Variant,
    VariantCollection,
)
from varcode.effects.effect_classes import HaplotypeEffect, Silent, Substitution


SHORT_SOURCE = "2024.06.11.bostongene.align.tcga.protocol.dr32Aligned.sorted"
LONG_SOURCE = "IPISRC044_T1_sclrs_ONT.tagged"
SHORT_NTF3 = SHORT_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_NTF3 = LONG_SOURCE + ".NTF3-chr12-5494381-compound"
LONG_EXOC4 = LONG_SOURCE + ".EXOC4-chr7-133274996"
MAP2_PREFIX = "isovar/osteosarc/figure_comparisons/corpus/context-evidence.json.gz#/map2/sources/"
SHORT_MAP2_T1 = MAP2_PREFIX + "11/original_records"
SHORT_MAP2_T2 = MAP2_PREFIX + "12/original_records"
LONG_MAP2_T2 = MAP2_PREFIX + "3/original_records"
MAP2_MEMBERS = (SHORT_MAP2_T1, SHORT_MAP2_T2, LONG_MAP2_T2)
MEMBERS = (SHORT_NTF3, LONG_NTF3, LONG_EXOC4) + MAP2_MEMBERS
MANIFEST_SHA256 = "b42dbce529cb5724153f035fca3783a41e673b31610a30851b2ef96d4ba43734"
REFERENCE_PATH = Path(__file__).parent / "data" / "sid_rna_effect_references.json"
REFERENCE_SHA256 = "d9e21cc96450d4a14e7b41d0900df48792c9f64b103b9ce01ab5cc841508d29e"

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


@pytest.fixture(scope="module")
def shared_phasing_reads(tmp_path_factory):
    osteosarc = pytest.importorskip("osteosarc", reason="install .[test-data]")
    pysam = pytest.importorskip("pysam", reason="install .[rna]")
    from osteosarc.records import read_records, read_template

    root = Path(osteosarc.fetch_bundle("openvax-v2"))
    manifest_bytes = (root / "manifest.json").read_bytes()
    assert hashlib.sha256(manifest_bytes).hexdigest() == MANIFEST_SHA256
    recipe = json.loads((root / "recipe.json").read_text())
    assert recipe["selection"]["supplementary"] is True
    compound = recipe["targets"]["NTF3-chr12-5494381-compound"]
    assert (compound["assembly"], compound["contig"], compound["position"],
            compound["ref"], compound["alt"]) == ("GRCh38", "chr12", 5494381, "AG", "GT")
    destination = tmp_path_factory.mktemp("sid-rna-phasing")
    paths = osteosarc.export_bundle("openvax-v2", destination, members=MEMBERS)
    manifest = json.loads(manifest_bytes)
    # The legacy MAP2 members are regional acquisitions. Include every other
    # original record of their selected templates that is available in the
    # bundle, including distant mates and supplementary placements.
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
    return paths, root, manifest


@pytest.mark.parametrize("member,records,templates,supplementary", [
    (SHORT_NTF3, 12, 6, 0),
    (LONG_NTF3, 15, 15, 0),
    (LONG_EXOC4, 40, 38, 2),
    (SHORT_MAP2_T1, 125, 78, 0),
    (SHORT_MAP2_T2, 26, 18, 0),
    (LONG_MAP2_T2, 33, 32, 1),
])
def test_shared_phasing_inputs_keep_complete_templates(
        shared_phasing_reads, member, records, templates, supplementary):
    from osteosarc.records import read_records

    paths, root, manifest = shared_phasing_reads
    exported = list(read_records(paths[member]))
    identities = Counter(record.digest for record in exported)
    member_identities = Counter(manifest["members"][member]["records"])
    assert not member_identities - identities
    if member not in MAP2_MEMBERS:
        assert identities == member_identities
    selected = {record.template for record in exported}
    source = manifest["members"][member]["source"]
    source_records = list(read_records(root / manifest["sources"][source]["bam"]))
    assert selected == {record.template for record in source_records
                        if record.digest in member_identities}
    # Compare to every record of each selected RG/QNAME in the bundle's source,
    # rather than refetching a locus and losing mates or distant split pieces.
    complete = Counter(record.digest for record in source_records if record.template in selected)
    assert identities == complete
    assert len(exported) == records
    assert len(selected) == templates
    assert sum(record.read.is_supplementary for record in exported) == supplementary
    if member == SHORT_NTF3:
        for template in selected:
            assert {record.read.flag & 0xc0 for record in exported
                    if record.template == template} == {0x40, 0x80}


@pytest.mark.parametrize("member,alleles,counts,support,phase", [
    (SHORT_NTF3, (NTF3, NTF3_COMPOUND_PARTNER), (3, 0, 0, 3), (3, 3), True),
    (LONG_NTF3, (NTF3, NTF3_COMPOUND_PARTNER), (9, 0, 0, 5), (9, 9), True),
    (LONG_NTF3, (NTF3, NTF3_COMPANION), (0, 9, 3, 2), (9, 3), False),
    (LONG_EXOC4, (EXOC4, EXOC4_COMPANION), (1, 0, 1, 0), (19, 3), None),
])
def test_real_sid_rna_fragment_evidence(
        shared_phasing_reads, member, alleles, counts, support, phase):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81)) for allele in alleles)
    source = RNAReadPhasingSource(
        str(paths[member]), variants=variants,
        # Q20/MAPQ20, five-base edge exclusion and SAM flag filters remain
        # at their defaults. Allow a 5% allele error rate for the ONT data.
        phasing_error_rate=0.01 if member == SHORT_NTF3 else 0.05)
    fetch = Mock(wraps=source._fetch_reads_for_variant)
    source._fetch_reads_for_variant = fetch
    try:
        assert source._phase_counts(*variants) == counts
        assert tuple(source.supports_variant(v) for v in variants) == support
        assert all(source.has_evidence(v) for v in variants)
        assert source.in_cis(*variants) is phase
        assert source.in_cis(*reversed(variants)) is phase
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(*variants) is phase
        for variant, partner in (variants, tuple(reversed(variants))):
            assert tuple(resolver.phased_partners(variant)) == ((partner,) if phase is True else ())
        assert tuple(source.supports_variant(v) for v in reversed(variants)) == support[::-1]
        assert fetch.call_count == len(variants)
    finally:
        source.close()


@pytest.fixture(scope="module")
def effect_references():
    """Offline RefSeq references and separately audited protein expectations."""
    payload = REFERENCE_PATH.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == REFERENCE_SHA256
    return json.loads(payload)


@pytest.mark.parametrize("gene", ["NTF3", "MAP2"])
def test_primary_effect_references_match_the_selected_transcripts(effect_references, gene):
    reference = effect_references["references"][gene]
    transcript = cached_release(81).transcript_by_id(reference["transcript_id"])
    assert transcript.coding_sequence == reference["cds"]
    assert transcript.protein_sequence == reference["protein"]
    assert hashlib.sha256(transcript.coding_sequence.encode()).hexdigest() == reference["cds_sha256"]
    assert hashlib.sha256(transcript.protein_sequence.encode()).hexdigest() == reference["protein_sha256"]
    window = reference["genomic_window"]
    # Independently acquired chromosome bases map to this exact CDS interval.
    start = min(transcript.start_codon_spliced_offsets)
    assert transcript.spliced_offset(window["start"]) - start == window["cds_offset"]
    assert transcript.sequence[transcript.spliced_offset(window["start"]):
                               transcript.spliced_offset(window["end"]) + 1] == window["sequence"]


def transcript_effects(effects, transcript_id, *, joint=False):
    return [effect for effect in effects
            if getattr(getattr(effect, "transcript", None), "id", None) == transcript_id
            and (len(getattr(effect, "variants", ())) > 1) is joint]


def individual_signatures(effects):
    """Snapshot all individual predictions, including other transcript models."""
    return sorted((str(effect.variant), getattr(getattr(effect, "transcript", None), "id", ""),
                   type(effect).__name__, effect.short_description, effect.mutant_protein_sequence)
                  for effect in effects if len(getattr(effect, "variants", ())) < 2)


def assert_primary_protein(effect, case, references):
    expected = references["expectations"][case]
    reference = references["references"][expected["gene"]]
    protein = (reference["protein"][:expected["protein_prefix_length"]]
               + expected["protein_replacement"]
               + reference["protein"][expected["protein_suffix_start"]:])
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
        cdnas = [outcome.mutant.cdna_sequence
                 for candidate in effect.realized_candidates for outcome in candidate.outcomes]
    assert cdnas
    start = min(effect.transcript.start_codon_spliced_offsets)
    delta = sum(len(alt) - len(ref) for _, ref, alt in expected["cds_edits"])
    for cdna in cdnas:
        cds = cdna[start:start + len(reference["cds"]) + delta]
        assert hashlib.sha256(cds.encode()).hexdigest() == expected["mutant_cds_sha256"]
        assert cds[expected["aa_pos"] * 3:expected["aa_pos"] * 3 + 3] == expected["alt_codon"]
        stop = expected["stop_codon_offset"]
        assert cds[stop:stop + 3] == expected["stop_codon"]


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("member", [SHORT_NTF3, LONG_NTF3])
def test_real_ntf3_cis_reads_predict_the_compound_codon(
        shared_phasing_reads, effect_references, annotator, member):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81))
                     for allele in (NTF3, NTF3_COMPOUND_PARTNER))
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = RNAReadPhasingSource(
        str(paths[member]), variants=variants,
        phasing_error_rate=0.01 if member == SHORT_NTF3 else 0.05)
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(*variants) is True
        effects = collection.effects(annotator=annotator, phase_resolver=resolver)
        assert individual_signatures(effects) == baseline
        transcript_id = effect_references["references"]["NTF3"]["transcript_id"]
        joint, = transcript_effects(effects, transcript_id, joint=True)
        assert isinstance(joint, Substitution if annotator == "transcript_model" else HaplotypeEffect)
        assert joint.variants == variants
        assert joint.phase_source == "rna_reads"
        assert joint.annotator == effects.annotator == annotator
        assert joint.annotator_version == effects.annotator_version
        assert_primary_protein(joint, "NTF3_compound", effect_references)
        assert_primary_cdna(joint, "NTF3_compound", effect_references)
        individuals = {effect.variant: effect for effect in transcript_effects(effects, transcript_id)}
        for variant, case in zip(variants, ("NTF3_first", "NTF3_second")):
            assert_primary_protein(individuals[variant], case, effect_references)
        restored = EffectCollection.from_json(effects.to_json())
        restored_joint, = transcript_effects(restored, transcript_id, joint=True)
        assert restored_joint.variants == variants
        assert restored_joint.phase_source == joint.phase_source
        assert_primary_protein(restored_joint, "NTF3_compound", effect_references)
    finally:
        source.close()


@pytest.mark.parametrize("annotator", ["fast", "protein_diff", "transcript_model"])
@pytest.mark.parametrize("alleles,cases,minimum,phase", [
    ((NTF3, NTF3_COMPANION), ("NTF3_first", "NTF3_companion"), 2, False),
    ((NTF3, NTF3_COMPOUND_PARTNER), ("NTF3_first", "NTF3_second"), 10, None),
])
def test_real_ntf3_trans_and_unknown_reads_keep_separate_predictions(
        shared_phasing_reads, effect_references, annotator, alleles, cases, minimum, phase):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81)) for allele in alleles)
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = RNAReadPhasingSource(
        str(paths[LONG_NTF3]), variants=variants, min_alt_reads=minimum, phasing_error_rate=0.05)
    try:
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(*variants) is phase  # Unknown must remain distinct from trans.
        assert all(not resolver.phased_partners(variant) for variant in variants)
        effects = collection.effects(annotator=annotator, phase_resolver=resolver)
        assert individual_signatures(effects) == baseline
        assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
        transcript_id = effect_references["references"]["NTF3"]["transcript_id"]
        individuals = {effect.variant: effect for effect in transcript_effects(effects, transcript_id)}
        assert set(individuals) == set(variants)
        for variant, case in zip(variants, cases):
            assert_primary_protein(individuals[variant], case, effect_references)
    finally:
        source.close()


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
@pytest.mark.parametrize("member,support,phase", [
    (SHORT_MAP2_T1, 4, True), (SHORT_MAP2_T2, 2, True), (LONG_MAP2_T2, 1, None),
])
def test_real_map2_anchored_reads_predict_the_complete_complex_allele(
        shared_phasing_reads, effect_references, annotator, member, support, phase):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81)) for allele in MAP2_PARTS)
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = RNAReadPhasingSource(
        str(paths[member]), variants=variants,
        phasing_error_rate=0.05 if member == LONG_MAP2_T2 else 0.01)
    try:
        source.register_haplotype(variants)
        assert tuple(source.supports_variant(variant) for variant in variants) == (support,) * 3
        resolver = MolecularPhaseResolver(source)
        for left, right in ((0, 1), (0, 2), (1, 2)):
            assert resolver.in_cis(variants[left], variants[right]) is phase
        effects = collection.effects(annotator=annotator, phase_resolver=resolver)
        assert individual_signatures(effects) == baseline
        transcript_id = effect_references["references"]["MAP2"]["transcript_id"]
        individuals = {effect.variant: effect for effect in transcript_effects(effects, transcript_id)}
        assert set(individuals) == set(variants)
        for variant, case in zip(variants, ("MAP2_first", "MAP2_second", "MAP2_deletion")):
            assert_primary_protein(individuals[variant], case, effect_references)
        if phase is None:
            assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
            assert not any(resolver.has_evidence(variant) for variant in variants)
        else:
            joint, = transcript_effects(effects, transcript_id, joint=True)
            assert joint.variants == variants
            assert joint.phase_source == "rna_reads"
            assert joint.annotator == effects.annotator == annotator
            assert_primary_protein(joint, "MAP2_compound", effect_references)
            assert_primary_cdna(joint, "MAP2_compound", effect_references)
            restored_joint, = transcript_effects(
                EffectCollection.from_json(effects.to_json()), transcript_id, joint=True)
            assert restored_joint.variants == variants
            assert restored_joint.phase_source == joint.phase_source
            assert_primary_protein(restored_joint, "MAP2_compound", effect_references)
            assert_primary_cdna(restored_joint, "MAP2_compound", effect_references)
    finally:
        source.close()


@pytest.mark.parametrize("annotator", ["fast", "protein_diff"])
@pytest.mark.parametrize("member,support", [
    (SHORT_MAP2_T1, 4), (SHORT_MAP2_T2, 2), (LONG_MAP2_T2, 1),
])
def test_real_map2_rna_does_not_conflate_competing_dna_alleles(
        shared_phasing_reads, effect_references, annotator, member, support):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81))
                     for allele in (MAP2_COMPOUND, MAP2_PARTS[-1], MAP2_HISTORICAL))
    collection = VariantCollection(variants)
    baseline = individual_signatures(collection.effects(annotator=annotator))
    source = RNAReadPhasingSource(
        str(paths[member]), variants=variants,
        phasing_error_rate=0.05 if member == LONG_MAP2_T2 else 0.01)
    try:
        for variant in variants:
            source.register_haplotype([variant])
        assert tuple(source.supports_variant(variant) for variant in variants) == (support, 0, 0)
        resolver = MolecularPhaseResolver(source)
        assert resolver.in_cis(variants[0], variants[1]) is None
        assert resolver.in_cis(variants[0], variants[2]) is None
        effects = collection.effects(annotator=annotator, phase_resolver=resolver)
        assert individual_signatures(effects) == baseline
        assert not any(len(getattr(effect, "variants", ())) > 1 for effect in effects)
        transcript_id = effect_references["references"]["MAP2"]["transcript_id"]
        individuals = {effect.variant: effect for effect in transcript_effects(effects, transcript_id)}
        assert set(individuals) == set(variants)
        for variant, case in zip(variants, ("MAP2_compound", "MAP2_deletion", "MAP2_historical")):
            assert_primary_protein(individuals[variant], case, effect_references)
    finally:
        source.close()


def test_real_sid_cis_call_requires_enough_fragments(shared_phasing_reads):
    paths, _, _ = shared_phasing_reads
    variants = tuple(Variant(*allele, genome=cached_release(81))
                     for allele in (NTF3, NTF3_COMPOUND_PARTNER))
    source = RNAReadPhasingSource(
        str(paths[LONG_NTF3]), variants=variants, min_alt_reads=10, phasing_error_rate=0.05)
    try:
        assert tuple(source.supports_variant(v) for v in variants) == (9, 9)
        assert not any(source.has_evidence(v) for v in variants)
        assert source.in_cis(*variants) is None
        assert MolecularPhaseResolver(source).in_cis(*variants) is None
    finally:
        source.close()


@pytest.mark.parametrize("skip_supplementary", [True, False])
def test_real_sid_supplementary_records_obey_the_evidence_filter(
        shared_phasing_reads, skip_supplementary):
    import pysam

    paths, _, _ = shared_phasing_reads
    variant = Variant(*EXOC4, genome=cached_release(81))
    with pysam.AlignmentFile(paths[LONG_EXOC4], "rb") as bam:
        supplementary = [read for read in bam if read.is_supplementary]
    assert len(supplementary) == 2
    source = RNAReadPhasingSource(
        str(paths[LONG_EXOC4]), skip_supplementary=skip_supplementary, phasing_error_rate=0.05)
    try:
        assert [source._read_allele(read, variant) for read in supplementary] == (
            [None, None] if skip_supplementary else ["ref", "ref"])
        assert source.supports_variant(variant) == 19
    finally:
        source.close()
