#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Predict variant consequences with explicit transcript and evidence context.

Start with :class:`Variant`, :func:`load_vcf`, and :meth:`Variant.effects`.
The default annotator is ``transcript_model``; ``fast`` is an explicit compatibility option.
Public imports below include compatibility names retained for existing callers.
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING, Any

from .version import __version__

if TYPE_CHECKING:
    from .annotators import (
        EffectAnnotator as EffectAnnotator,
        FastEffectAnnotator as FastEffectAnnotator,
        ProteinDiffEffectAnnotator as ProteinDiffEffectAnnotator,
        RealizedEffectAnnotator as RealizedEffectAnnotator,
        TranscriptModelEffectAnnotator as TranscriptModelEffectAnnotator,
        get_annotator as get_annotator,
        get_default_annotator as get_default_annotator,
        register_annotator as register_annotator,
        registered_annotators as registered_annotators,
        resolve_annotator as resolve_annotator,
        set_default_annotator as set_default_annotator,
        use_annotator as use_annotator,
    )
    from .contigs import (
        ContigResolution as ContigResolution,
        ContigResolver as ContigResolver,
        normalize_contig_name as normalize_contig_name,
        resolve_contig_name as resolve_contig_name,
    )
    from .cryptic_exons import (
        enumerate_candidates as enumerate_cryptic_exon_candidates,
        enumerate_from_structural_variant as enumerate_from_structural_variant,
        score_acceptor as score_acceptor,
        score_donor as score_donor,
    )
    from .effect_candidates import (
        EffectCandidate as EffectCandidate,
        candidates_from_effects as candidates_from_effects,
    )
    from .effect_traversal import (
        Candidate as Candidate,
        CandidatePath as CandidatePath,
        get_effect_candidates as get_effect_candidates,
        iter_candidate_effects as iter_candidate_effects,
        iter_candidate_paths as iter_candidate_paths,
    )
    from .effects import (
        EffectCollection as EffectCollection,
        EffectReference as EffectReference,
        GermlineAlleleOverlap as GermlineAlleleOverlap,
        HypothesisLimit as HypothesisLimit,
        MultiOutcomeEffect as MultiOutcomeEffect,
        MutationEffect as MutationEffect,
        NonsilentCodingMutation as NonsilentCodingMutation,
        Unresolved as Unresolved,
        classification_priority as classification_priority,
        effect_collection_from_json as effect_collection_from_json,
        effect_collection_from_serializable_repr as effect_collection_from_serializable_repr,
        effect_collection_to_json as effect_collection_to_json,
        effect_from_json as effect_from_json,
        effect_from_serializable_repr as effect_from_serializable_repr,
        effect_priority as effect_priority,
        effect_serialization_classes as effect_serialization_classes,
        effect_to_json as effect_to_json,
        top_priority_effect as top_priority_effect,
    )
    from .effects.effect_classes import (
        CrypticAcceptor as CrypticAcceptor,
        CrypticDonor as CrypticDonor,
        ExonSkipping as ExonSkipping,
        ExonicSpliceSite as ExonicSpliceSite,
        IntronRetention as IntronRetention,
        IntronicSpliceSite as IntronicSpliceSite,
        NormalSplicing as NormalSplicing,
        SpliceAcceptor as SpliceAcceptor,
        SpliceDonor as SpliceDonor,
        SpliceMechanismEffect as SpliceMechanismEffect,
        SpliceSite as SpliceSite,
    )
    from .errors import (
        HypothesisLimitError as HypothesisLimitError,
        InvalidContigError as InvalidContigError,
        ReferenceMismatchError as ReferenceMismatchError,
        SampleNotFoundError as SampleNotFoundError,
    )
    from .exacto import load_exacto_fusions as load_exacto_fusions
    from .genome import Genome as Genome
    from .genome_sequence import (
        reference_base as reference_base,
        reference_range as reference_range,
    )
    from .genotype import (
        Genotype as Genotype,
        Zygosity as Zygosity,
    )
    from .germline import (
        Completeness as Completeness,
        GenomeBuildMismatchError as GenomeBuildMismatchError,
        GermlineContext as GermlineContext,
        PhaseHypothesis as PhaseHypothesis,
        PhasePartition as PhasePartition,
        apply_germline_to_transcript as apply_germline_to_transcript,
        default_germline_window as default_germline_window,
        detect_germline_overlap as detect_germline_overlap,
        detect_loh as detect_loh,
        enumerate_phase_hypotheses as enumerate_phase_hypotheses,
        partition_germline_by_phase as partition_germline_by_phase,
        predict_germline_aware_effect as predict_germline_aware_effect,
    )
    from .maf import (
        load_maf as load_maf,
        load_maf_dataframe as load_maf_dataframe,
    )
    from .mutant_transcript import (
        AssembledAllele as AssembledAllele,
        MutantTranscript as MutantTranscript,
        ReferenceSegment as ReferenceSegment,
        TranscriptEdit as TranscriptEdit,
        apply_variant_to_transcript as apply_variant_to_transcript,
        apply_variants_to_transcript as apply_variants_to_transcript,
    )
    from .orf import (
        ORFAnnotation as ORFAnnotation,
        ORFOutcome as ORFOutcome,
        ORFProduct as ORFProduct,
        ORFProductAnnotation as ORFProductAnnotation,
        annotate_orf_product as annotate_orf_product,
        annotate_orf_rna_product as annotate_orf_rna_product,
        annotate_orf_variants as annotate_orf_variants,
        load_orf_catalog as load_orf_catalog,
        predict_orf_effects as predict_orf_effects,
    )
    from .orf_discovery import (
        DiscoveredORF as DiscoveredORF,
        ORFDiscovery as ORFDiscovery,
        ORFDiscoveryAssessment as ORFDiscoveryAssessment,
        ORFDiscoveryPolicy as ORFDiscoveryPolicy,
        discover_created_orfs as discover_created_orfs,
        discover_created_orfs_in_products as discover_created_orfs_in_products,
    )
    from .phasing import (
        MolecularPhaseResolver as MolecularPhaseResolver,
        MutantTranscriptSource as MutantTranscriptSource,
        ReadPhaseResolver as ReadPhaseResolver,
        ReadPhasingSource as ReadPhasingSource,
        VCFPhaseResolver as VCFPhaseResolver,
        apply_phase_resolver_to_effects as apply_phase_resolver_to_effects,
        query_in_cis as query_in_cis,
    )
    from .proteome import to_proteome as to_proteome
    from .realized_effects import predict_realized_effect as predict_realized_effect
    from .reference_completion import (
        reference_completion_hypotheses as reference_completion_hypotheses,
    )
    from .rna_evidence import (
        NullRNAEvidenceResolver as NullRNAEvidenceResolver,
        RNAEvidence as RNAEvidence,
        RNAEvidenceResolver as RNAEvidenceResolver,
        apply_rna_evidence_to_effects as apply_rna_evidence_to_effects,
        make_fusion_outcome as make_fusion_outcome,
        make_rna_outcome as make_rna_outcome,
    )
    from .rna_read_phasing import RNAReadPhasingSource as RNAReadPhasingSource
    from .sample_identity import (
        SampleCheckConfig as SampleCheckConfig,
        SampleFingerprint as SampleFingerprint,
        SampleSpec as SampleSpec,
        check_sample_identity as check_sample_identity,
        compare_samples as compare_samples,
        load_vcf_samples as load_vcf_samples,
    )
    from .splice_outcomes import (
        SpliceCandidateRNAEvidence as SpliceCandidateRNAEvidence,
        SpliceOutcomeSet as SpliceOutcomeSet,
    )
    from .structural_variant import (
        SV_TYPES as SV_TYPES,
        StructuralVariant as StructuralVariant,
    )
    from .sv_allele_parser import parse_symbolic_alt as parse_symbolic_alt
    from .transcript_layout import CrypticSplicePolicy as CrypticSplicePolicy
    from .transcript_model import (
        predict_transcript_model_effect as predict_transcript_model_effect,
    )
    from .transforms import (
        combine_cis_snvs as combine_cis_snvs,
    )
    from .uorf import (
        ORFTranscript as ORFTranscript,
        UORFAnnotation as UORFAnnotation,
        annotate_uorf_product as annotate_uorf_product,
        annotate_uorf_variants as annotate_uorf_variants,
        load_uorf_catalog as load_uorf_catalog,
        predict_uorf_effects as predict_uorf_effects,
    )
    from .variant import Variant as Variant
    from .variant_collection import (
        SpanningDeletionRecord as SpanningDeletionRecord,
        VariantCollection as VariantCollection,
    )
    from .vcf import (
        load_vcf as load_vcf,
        load_vcf_fast as load_vcf_fast,
    )


_LAZY_EXPORTS: dict[str, tuple[str, str]] = {
    # annotators
    "EffectAnnotator": (".annotators", "EffectAnnotator"),
    "FastEffectAnnotator": (".annotators", "FastEffectAnnotator"),
    "ProteinDiffEffectAnnotator": (".annotators", "ProteinDiffEffectAnnotator"),
    "RealizedEffectAnnotator": (".annotators", "RealizedEffectAnnotator"),
    "TranscriptModelEffectAnnotator": (".annotators", "TranscriptModelEffectAnnotator"),
    "get_annotator": (".annotators", "get_annotator"),
    "get_default_annotator": (".annotators", "get_default_annotator"),
    "register_annotator": (".annotators", "register_annotator"),
    "registered_annotators": (".annotators", "registered_annotators"),
    "resolve_annotator": (".annotators", "resolve_annotator"),
    "set_default_annotator": (".annotators", "set_default_annotator"),
    "use_annotator": (".annotators", "use_annotator"),
    # exceptions
    "InvalidContigError": (".errors", "InvalidContigError"),
    "HypothesisLimitError": (".errors", "HypothesisLimitError"),
    "ReferenceMismatchError": (".errors", "ReferenceMismatchError"),
    "SampleNotFoundError": (".errors", "SampleNotFoundError"),
    # genome / contigs / reference DNA
    "Genome": (".genome", "Genome"),
    "ContigResolution": (".contigs", "ContigResolution"),
    "ContigResolver": (".contigs", "ContigResolver"),
    "normalize_contig_name": (".contigs", "normalize_contig_name"),
    "resolve_contig_name": (".contigs", "resolve_contig_name"),
    "reference_base": (".genome_sequence", "reference_base"),
    "reference_range": (".genome_sequence", "reference_range"),
    # germline
    "GenomeBuildMismatchError": (".germline", "GenomeBuildMismatchError"),
    "GermlineContext": (".germline", "GermlineContext"),
    "Completeness": (".germline", "Completeness"),
    "PhaseHypothesis": (".germline", "PhaseHypothesis"),
    "PhasePartition": (".germline", "PhasePartition"),
    "apply_germline_to_transcript": (".germline", "apply_germline_to_transcript"),
    "default_germline_window": (".germline", "default_germline_window"),
    "detect_loh": (".germline", "detect_loh"),
    "detect_germline_overlap": (".germline", "detect_germline_overlap"),
    "enumerate_phase_hypotheses": (".germline", "enumerate_phase_hypotheses"),
    "partition_germline_by_phase": (".germline", "partition_germline_by_phase"),
    "predict_germline_aware_effect": (".germline", "predict_germline_aware_effect"),
    # genotype
    "Genotype": (".genotype", "Genotype"),
    "Zygosity": (".genotype", "Zygosity"),
    # mutant_transcript
    "AssembledAllele": (".mutant_transcript", "AssembledAllele"),
    "MutantTranscript": (".mutant_transcript", "MutantTranscript"),
    "ReferenceSegment": (".mutant_transcript", "ReferenceSegment"),
    "TranscriptEdit": (".mutant_transcript", "TranscriptEdit"),
    "apply_variant_to_transcript": (".mutant_transcript", "apply_variant_to_transcript"),
    "apply_variants_to_transcript": (".mutant_transcript", "apply_variants_to_transcript"),
    # cryptic_exons
    "enumerate_cryptic_exon_candidates": (".cryptic_exons", "enumerate_candidates"),
    "enumerate_from_structural_variant": (".cryptic_exons", "enumerate_from_structural_variant"),
    "score_acceptor": (".cryptic_exons", "score_acceptor"),
    "score_donor": (".cryptic_exons", "score_donor"),
    # effect candidates & traversal
    "EffectCandidate": (".effect_candidates", "EffectCandidate"),
    "candidates_from_effects": (".effect_candidates", "candidates_from_effects"),
    "Candidate": (".effect_traversal", "Candidate"),
    "CandidatePath": (".effect_traversal", "CandidatePath"),
    "get_effect_candidates": (".effect_traversal", "get_effect_candidates"),
    "iter_candidate_paths": (".effect_traversal", "iter_candidate_paths"),
    "iter_candidate_effects": (".effect_traversal", "iter_candidate_effects"),
    # phasing & RNA evidence
    "MolecularPhaseResolver": (".phasing", "MolecularPhaseResolver"),
    "MutantTranscriptSource": (".phasing", "MutantTranscriptSource"),
    "ReadPhaseResolver": (".phasing", "ReadPhaseResolver"),
    "ReadPhasingSource": (".phasing", "ReadPhasingSource"),
    "VCFPhaseResolver": (".phasing", "VCFPhaseResolver"),
    "apply_phase_resolver_to_effects": (".phasing", "apply_phase_resolver_to_effects"),
    "query_in_cis": (".phasing", "query_in_cis"),
    "RNAReadPhasingSource": (".rna_read_phasing", "RNAReadPhasingSource"),
    "NullRNAEvidenceResolver": (".rna_evidence", "NullRNAEvidenceResolver"),
    "RNAEvidence": (".rna_evidence", "RNAEvidence"),
    "RNAEvidenceResolver": (".rna_evidence", "RNAEvidenceResolver"),
    "apply_rna_evidence_to_effects": (".rna_evidence", "apply_rna_evidence_to_effects"),
    "make_fusion_outcome": (".rna_evidence", "make_fusion_outcome"),
    "make_rna_outcome": (".rna_evidence", "make_rna_outcome"),
    "load_exacto_fusions": (".exacto", "load_exacto_fusions"),
    "reference_completion_hypotheses": (".reference_completion", "reference_completion_hypotheses"),
    "predict_realized_effect": (".realized_effects", "predict_realized_effect"),
    "CrypticSplicePolicy": (".transcript_layout", "CrypticSplicePolicy"),
    "predict_transcript_model_effect": (".transcript_model", "predict_transcript_model_effect"),
    # splice outcomes & splice effect classes
    "SpliceCandidateRNAEvidence": (".splice_outcomes", "SpliceCandidateRNAEvidence"),
    "SpliceOutcomeSet": (".splice_outcomes", "SpliceOutcomeSet"),
    "CrypticAcceptor": (".effects.effect_classes", "CrypticAcceptor"),
    "CrypticDonor": (".effects.effect_classes", "CrypticDonor"),
    "ExonSkipping": (".effects.effect_classes", "ExonSkipping"),
    "ExonicSpliceSite": (".effects.effect_classes", "ExonicSpliceSite"),
    "IntronRetention": (".effects.effect_classes", "IntronRetention"),
    "IntronicSpliceSite": (".effects.effect_classes", "IntronicSpliceSite"),
    "NormalSplicing": (".effects.effect_classes", "NormalSplicing"),
    "SpliceAcceptor": (".effects.effect_classes", "SpliceAcceptor"),
    "SpliceDonor": (".effects.effect_classes", "SpliceDonor"),
    "SpliceMechanismEffect": (".effects.effect_classes", "SpliceMechanismEffect"),
    "SpliceSite": (".effects.effect_classes", "SpliceSite"),
    # structural variants & variants
    "SV_TYPES": (".structural_variant", "SV_TYPES"),
    "StructuralVariant": (".structural_variant", "StructuralVariant"),
    "parse_symbolic_alt": (".sv_allele_parser", "parse_symbolic_alt"),
    "Variant": (".variant", "Variant"),
    "SpanningDeletionRecord": (".variant_collection", "SpanningDeletionRecord"),
    "VariantCollection": (".variant_collection", "VariantCollection"),
    "combine_cis_snvs": (".transforms", "combine_cis_snvs"),
    # file loading
    "load_maf": (".maf", "load_maf"),
    "load_maf_dataframe": (".maf", "load_maf_dataframe"),
    "load_vcf": (".vcf", "load_vcf"),
    "load_vcf_fast": (".vcf", "load_vcf_fast"),
    # effects & serialization
    "effect_priority": (".effects", "effect_priority"),
    "classification_priority": (".effects", "classification_priority"),
    "top_priority_effect": (".effects", "top_priority_effect"),
    "EffectCollection": (".effects", "EffectCollection"),
    "EffectReference": (".effects", "EffectReference"),
    "effect_collection_from_json": (".effects", "effect_collection_from_json"),
    "effect_collection_from_serializable_repr": (".effects", "effect_collection_from_serializable_repr"),
    "effect_collection_to_json": (".effects", "effect_collection_to_json"),
    "effect_from_json": (".effects", "effect_from_json"),
    "effect_from_serializable_repr": (".effects", "effect_from_serializable_repr"),
    "effect_serialization_classes": (".effects", "effect_serialization_classes"),
    "effect_to_json": (".effects", "effect_to_json"),
    "MultiOutcomeEffect": (".effects", "MultiOutcomeEffect"),
    "MutationEffect": (".effects", "MutationEffect"),
    "NonsilentCodingMutation": (".effects", "NonsilentCodingMutation"),
    "HypothesisLimit": (".effects", "HypothesisLimit"),
    "Unresolved": (".effects", "Unresolved"),
    "GermlineAlleleOverlap": (".effects", "GermlineAlleleOverlap"),
    # proteome & sample identity
    "to_proteome": (".proteome", "to_proteome"),
    "SampleCheckConfig": (".sample_identity", "SampleCheckConfig"),
    "SampleSpec": (".sample_identity", "SampleSpec"),
    "SampleFingerprint": (".sample_identity", "SampleFingerprint"),
    "load_vcf_samples": (".sample_identity", "load_vcf_samples"),
    "compare_samples": (".sample_identity", "compare_samples"),
    "check_sample_identity": (".sample_identity", "check_sample_identity"),
    # ORF / uORF / discovery
    "ORFTranscript": (".uorf", "ORFTranscript"),
    "UORFAnnotation": (".uorf", "UORFAnnotation"),
    "annotate_uorf_variants": (".uorf", "annotate_uorf_variants"),
    "annotate_uorf_product": (".uorf", "annotate_uorf_product"),
    "predict_uorf_effects": (".uorf", "predict_uorf_effects"),
    "load_uorf_catalog": (".uorf", "load_uorf_catalog"),
    "ORFAnnotation": (".orf", "ORFAnnotation"),
    "ORFOutcome": (".orf", "ORFOutcome"),
    "ORFProduct": (".orf", "ORFProduct"),
    "ORFProductAnnotation": (".orf", "ORFProductAnnotation"),
    "annotate_orf_variants": (".orf", "annotate_orf_variants"),
    "annotate_orf_rna_product": (".orf", "annotate_orf_rna_product"),
    "annotate_orf_product": (".orf", "annotate_orf_product"),
    "predict_orf_effects": (".orf", "predict_orf_effects"),
    "load_orf_catalog": (".orf", "load_orf_catalog"),
    "DiscoveredORF": (".orf_discovery", "DiscoveredORF"),
    "ORFDiscovery": (".orf_discovery", "ORFDiscovery"),
    "ORFDiscoveryAssessment": (".orf_discovery", "ORFDiscoveryAssessment"),
    "ORFDiscoveryPolicy": (".orf_discovery", "ORFDiscoveryPolicy"),
    "discover_created_orfs": (".orf_discovery", "discover_created_orfs"),
    "discover_created_orfs_in_products": (".orf_discovery", "discover_created_orfs_in_products"),
}

_LAZY_SUBMODULES: frozenset[str] = frozenset({
    "annotators",
    "contigs",
    "cryptic_exons",
    "effect_candidates",
    "effect_traversal",
    "effects",
    "errors",
    "exacto",
    "extra_orfs",
    "genome",
    "genome_sequence",
    "genotype",
    "germline",
    "maf",
    "mutant_transcript",
    "orf",
    "orf_discovery",
    "phasing",
    "proteome",
    "realized_effects",
    "reference",
    "reference_completion",
    "rna_evidence",
    "rna_read_phasing",
    "sample_identity",
    "splice_outcomes",
    "structural_variant",
    "sv_allele_parser",
    "transcript_model",
    "transforms",
    "uorf",
    "variant",
    "variant_collection",
    "vcf",
})


__all__ = [
    "ORFAnnotation", "ORFOutcome", "ORFProduct", "ORFProductAnnotation",
    "annotate_orf_variants", "annotate_orf_rna_product", "annotate_orf_product",
    "predict_orf_effects", "load_orf_catalog",
    "ORFTranscript", "UORFAnnotation", "annotate_uorf_variants",
    "annotate_uorf_product", "predict_uorf_effects", "load_uorf_catalog",
    "SampleCheckConfig", "SampleSpec", "SampleFingerprint",
    "load_vcf_samples", "compare_samples", "check_sample_identity",
    "__version__",
    "to_proteome",
    "Candidate", "CandidatePath", "get_effect_candidates",
    "iter_candidate_paths", "iter_candidate_effects",

    # basic classes
    "Variant",
    "EffectCollection",
    "SpanningDeletionRecord",
    "VariantCollection",
    "combine_cis_snvs",

    # genotype / zygosity
    "Genotype",
    "Zygosity",

    # splice-signal disruption effects — "where" the signal was hit
    # (all share the SpliceSite base)
    "SpliceSite",
    "SpliceDonor",
    "SpliceAcceptor",
    "IntronicSpliceSite",
    "ExonicSpliceSite",

    # splice outcome possibility set (openvax/varcode#262, #382)
    # — "what" the spliceosome does (the protein-affecting effects)
    "SpliceOutcomeSet",
    "SpliceCandidateRNAEvidence",
    "SpliceMechanismEffect",
    "NormalSplicing",
    "ExonSkipping",
    "IntronRetention",
    "CrypticDonor",
    "CrypticAcceptor",

    # MutantTranscript data model + pluggable annotators (openvax/varcode#271)
    "AssembledAllele",
    "MutantTranscript",
    "ReferenceSegment",
    "TranscriptEdit",
    "apply_variant_to_transcript",
    "apply_variants_to_transcript",

    # Unified multi-outcome type (openvax/varcode#299)
    "EffectCandidate",
    "candidates_from_effects",
    "CrypticSplicePolicy",
    "predict_realized_effect",
    "predict_transcript_model_effect",

    # Structural variants (openvax/varcode#252 / #264)
    "StructuralVariant",
    "SV_TYPES",
    "parse_symbolic_alt",

    # Cryptic-exon candidate enumerator (pluggable scorer hook)
    "enumerate_cryptic_exon_candidates",
    "enumerate_from_structural_variant",
    "score_donor",
    "score_acceptor",

    # Phase resolvers (openvax/varcode#269)
    "MolecularPhaseResolver",
    "MutantTranscriptSource",
    "ReadPhaseResolver",
    "ReadPhasingSource",
    "RNAReadPhasingSource",
    "VCFPhaseResolver",
    "apply_phase_resolver_to_effects",
    "query_in_cis",

    # RNA-evidence resolver (openvax/varcode#259)
    "RNAEvidenceResolver",
    "RNAEvidence",
    "NullRNAEvidenceResolver",
    "apply_rna_evidence_to_effects",
    "make_rna_outcome",
    "make_fusion_outcome",
    "load_exacto_fusions",
    "reference_completion_hypotheses",

    "EffectAnnotator",
    "FastEffectAnnotator",
    "ProteinDiffEffectAnnotator",
    "RealizedEffectAnnotator",
    "TranscriptModelEffectAnnotator",
    "get_annotator",
    "get_default_annotator",
    "register_annotator",
    "registered_annotators",
    "resolve_annotator",
    "set_default_annotator",
    "use_annotator",

    # effects
    "EffectReference",
    "effect_collection_from_json",
    "effect_collection_from_serializable_repr",
    "effect_collection_to_json",
    "effect_from_json",
    "effect_from_serializable_repr",
    "effect_serialization_classes",
    "effect_to_json",
    "effect_priority",
    "classification_priority",
    "top_priority_effect",
    "MultiOutcomeEffect",
    "MutationEffect",
    "HypothesisLimit",
    "Unresolved",
    "NonsilentCodingMutation",

    # exceptions
    "InvalidContigError",
    "ReferenceMismatchError",
    "SampleNotFoundError",
    "GenomeBuildMismatchError",
    "HypothesisLimitError",

    # Germline-aware annotation (openvax/varcode#268)
    "GermlineContext",
    "Completeness",
    "PhaseHypothesis",
    "PhasePartition",
    "apply_germline_to_transcript",
    "default_germline_window",
    "detect_loh",
    "detect_germline_overlap",
    "GermlineAlleleOverlap",
    "enumerate_phase_hypotheses",
    "partition_germline_by_phase",
    "predict_germline_aware_effect",

    # Genome wrapper + tiered reference lookup (openvax/varcode#372)
    "Genome",
    "reference_base",
    "reference_range",

    # Reference-vocabulary name resolution (no coordinate conversion)
    "ContigResolution",
    "ContigResolver",
    "normalize_contig_name",
    "resolve_contig_name",

    # file loading
    "load_maf",
    "load_maf_dataframe",
    "load_vcf",
    "load_vcf_fast",

    # created-ORF discovery
    "DiscoveredORF",
    "ORFDiscovery",
    "ORFDiscoveryAssessment",
    "ORFDiscoveryPolicy",
    "discover_created_orfs",
    "discover_created_orfs_in_products",
]


def __getattr__(name: str) -> Any:
    target = _LAZY_EXPORTS.get(name)
    if target is not None:
        module_name, attr_name = target
        module = importlib.import_module(module_name, __name__)
        value = getattr(module, attr_name)
        globals()[name] = value
        return value
    if name in _LAZY_SUBMODULES:
        module = importlib.import_module(f".{name}", __name__)
        globals()[name] = module
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
