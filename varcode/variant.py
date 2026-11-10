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

from __future__ import annotations

from typing import TYPE_CHECKING

from serializable import Serializable

if TYPE_CHECKING:
    from pyensembl import Transcript

    from .annotators import EffectAnnotator
    from .effects import EffectCollection, MutationEffect

from . import circular_origin
from .contigs import (
    _HG19_MITOCHONDRION_NOTE, _annotation_contig_name, _convert_ucsc_contig_name,
    _is_hg19_mitochondrion)
from .nucleotides import (
    normalize_nucleotide_string,
    STANDARD_NUCLEOTIDES,
    is_purine
)
from .pinned_references import (
    is_pinned_assembly, mitochondrion_note, names_rcrs, pinned_contig_name)
from .reference import (
    _split_release, ensembl_to_ucsc_reference_names, infer_genome, infer_reference_name)
from .string_helpers import trim_shared_flanking_strings
from .errors import InvalidContigError
from .effects import (
    predict_variant_effects,
    predict_variant_effect_on_transcript
)

class Variant(Serializable):
    contig: str

    __slots__ = (
        "contig",
        "start",
        "end",
        "ref",
        "alt",
        "original_genome",
        "original_reference_name",
        "genome",
        "normalize_contig_names",
        "original_genome_was_ucsc",
        "convert_ucsc_contig_names",
        "chrm_is_rcrs",
        "allow_extended_nucleotides",
        "original_contig",
        "original_ref",
        "original_alt",
        "original_start",
        "_transcripts",
        "_genes",
        "_gene_ids",
        "_gene_names",
    )

    # any keywords which are renamed or removed should be added to this
    # dictionary to preserve the ability to deserialize old representations
    # of Variant objects
    _SERIALIZABLE_KEYWORD_ALIASES = {
        "normalize_contig_name": "normalize_contig_names"
    }
    
    def __init__(
            self,
            contig,
            start,
            ref,
            alt,
            genome=None,
            ensembl=None,
            allow_extended_nucleotides=False,
            normalize_contig_names=True,
            convert_ucsc_contig_names=None,
            *,
            original_reference_name=None,
            original_genome_was_ucsc=None,
            original_genome=None,
            chrm_is_rcrs=False):
        """
        Construct a Variant object.

        Parameters
        ----------
        contig : str
            Chromosome that this variant is on

        start : int
            1-based position on the chromosome of first reference nucleotide

        ref : str
            Reference nucleotide(s)

        alt : str
            Alternate nucleotide(s)

        genome : Genome, EnsemblRelease, or str, or int
            Name of reference genome, Ensembl release number, or object
            derived from pyensembl.Genome. Default to latest available release
            of GRCh38

        ensembl : Genome, EnsemblRelease, or str, or int (DEPRECATED)
            Previous name used instead of 'genome', the two arguments should
            be mutually exclusive.

        allow_extended_nucleotides : bool
            Extended nucleotides include 'Y' for pyrimidies or 'N' for any base

        normalize_contig_names : bool
            By default the contig name will be normalized by converting integers
            to strings (e.g. 1 -> "1"), and converting any letters after "chr"
            to uppercase (e.g. "chrx" -> "chrX"). If you don't want
            this behavior then pass normalize_contig_name=False.

        convert_ucsc_contig_names : bool, optional
            Setting this argument to True causes UCSC chromosome names to be
            coverted, such as "chr1" to "1". If the default value (None) is used
            then it defaults to whether or not a UCSC genome was pass in for
            the 'genome' argument. Under GRCh37, "chrM" and "M" become "chrM",
            not "MT": hg19's chrM is not the rCRS that Ensembl names MT, so
            such a variant can't be annotated, unless ``chrm_is_rcrs`` is set.
            A name of a genome ORFdata pins (e.g. "hs1" or
            "chm13v2.0_maskedY_rCRS.fa") converts too, by ORFdata's contig
            table: T2T-CHM13v2.0's "chr1", "CP068277.2" and "NC_060925.1" are
            "1".

        chrm_is_rcrs : bool, default False
            This input's "chrM" or "M" is the rCRS. Under GRCh37 it is then
            the rCRS that Ensembl names MT (as in GENCODE's GRCh37), not
            hg19's NC_001807, so UCSC conversion names it "MT". :func:`load_vcf`
            sets it when the header declares chrM at the rCRS length, 16,569
            bp. Under T2T-CHM13v2.0 it is the rCRS of the
            ``chm13v2.0_maskedY_rCRS`` analysis set, named "NC_012920.1" and
            annotated with Ensembl GRCh38's MT genes; without it, "chrM" is
            CHM13's own mitochondrion, "MT", which has no annotation. It is
            set when ``genome`` names such a file. No effect under other
            builds, whose chrM is the rCRS.

        original_reference_name, original_genome_was_ucsc, original_genome : optional
            Saved input provenance when reconstructing a pinned genome. These
            labels do not resolve another annotation or convert coordinates.
        """

        # first initialize the fields we use to cache lists of overlapping
        # pyensembl Gene and Transcript objects, or their properties such
        # as names/IDs
        self._genes = None
        self._transcripts = None
        self._gene_ids = None
        self._gene_names = None

        # Provenance: when this variant was produced by a
        # ``varcode.transforms`` operation (e.g. ``pair_breakends``,
        # ``combine_cis_snvs``), the source rows are recorded here.
        # Empty tuple for variants loaded directly from a VCF/MAF.
        # Not part of hash/equality — identity is still
        # (contig, start, ref, alt, reference_name).
        self.source_variants: tuple = ()

        # store the options which affect how properties of this variant
        # may be changed/transformed
        self.normalize_contig_names = normalize_contig_names
        self.allow_extended_nucleotides = allow_extended_nucleotides

        # if genome not specified, try the old name 'ensembl'
        # if ensembl is also None, then default to "GRCh38"
        if genome is None and ensembl is None:
            genome = "GRCh38"
        elif genome is None:
            genome = ensembl


        # user might supply Ensembl release as an integer, reference name,
        # or pyensembl.Genome object
        self.original_genome = genome
        self.genome, self.original_genome_was_ucsc = infer_genome(genome)

        self.reference_name = self.genome.reference_name
        # A file of a genome ORFdata pins, e.g. "chm13v2.0_maskedY_rCRS.fa".
        pinned = (infer_reference_name(_split_release(genome)[0])
                  if self.original_genome_was_ucsc and isinstance(genome, str)
                  and is_pinned_assembly(self.reference_name) else None)
        if pinned is not None:
            self.original_reference_name = pinned
        elif self.original_genome_was_ucsc:
            # Historical dated datasets use accessions as reference names.
            self.original_reference_name = ensembl_to_ucsc_reference_names[
                getattr(self.genome, "assembly", self.reference_name)]
        else:
            self.original_reference_name = self.reference_name

        if original_reference_name is not None:
            self.original_reference_name = original_reference_name
        if original_genome_was_ucsc is not None:
            self.original_genome_was_ucsc = original_genome_was_ucsc
        if original_genome is not None:
            self.original_genome = original_genome

        self.original_contig = contig

        if convert_ucsc_contig_names is None:
            self.convert_ucsc_contig_names = self.original_genome_was_ucsc
        else:
            self.convert_ucsc_contig_names = convert_ucsc_contig_names
        # The rCRS analysis set of T2T-CHM13 carries the rCRS as chrM.
        self.chrm_is_rcrs = bool(chrm_is_rcrs) or (pinned is not None and names_rcrs(pinned))

        self.contig = self._normalize_contig_name(contig)

        if ref != alt and ref in STANDARD_NUCLEOTIDES and alt in STANDARD_NUCLEOTIDES:
            # Optimization for common case.
            self.original_ref = self.ref = ref
            self.original_alt = self.alt = alt
            self.original_start = self.start = self.end = int(start)
            return

        # we want to preserve the ref/alt/pos both as they appeared in the
        # original VCF or MAF file but also normalize variants to get rid
        # of shared prefixes/suffixes between the ref and alt nucleotide
        # strings e.g. g.10 CTT>T can be normalized into g.10delCT
        #
        # The normalized variant properties go into fields
        #    Variant.{original_ref, original_alt, original_pos}
        # whereas the trimmed fields are:
        #    Variant.{ref, alt, start, end}

        # the original entries must preserve the number of nucleotides in
        # ref and alt but we still want to normalize e.g. '-' and '.' into ''
        self.original_ref = normalize_nucleotide_string(
            ref,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_alt = normalize_nucleotide_string(
            alt,
            allow_extended_nucleotides=allow_extended_nucleotides)
        self.original_start = int(start)

        # normalize the variant by trimming any shared prefix or suffix
        # between ref and alt nucleotide sequences and then
        # offset the variant position in a strand-dependent manner
        (trimmed_ref, trimmed_alt, prefix, _) = \
            trim_shared_flanking_strings(self.original_ref, self.original_alt)

        self.ref = trimmed_ref
        self.alt = trimmed_alt

        if len(trimmed_ref) == 0:
            # insertions must be treated differently since the meaning of a
            # position for an insertion is:
            #   "insert the alt nucleotides after this position"
            #
            # Aside: what if both trimmed ref and alt strings are empty?
            # This means we had a "null" variant, probably from a VCF
            # generated by force-calling mutations which weren't actually
            # found in the sample.
            # Null variants are interepted as inserting zero nucleotides
            # after the whole reference sequence.
            #
            # Start and end both are base-1 nucleotide position before
            # insertion.
            self.start = self.original_start + max(0, len(prefix) - 1)
            self.end = self.start
        else:
            # for substitutions and deletions the [start:end] interval is
            # an inclusive selection of reference nucleotides
            self.start = self.original_start + len(prefix)
            self.end = self.start + len(trimmed_ref) - 1


    @property
    def ensembl(self):
        """
        Deprecated alias for Variant.genome

        Returns
        -------
        pyensembl.Genome
        """
        return self.genome

    def __str__(self):
        return "Variant(contig='%s', start=%d, ref='%s', alt='%s', reference_name='%s')" % (
            self.contig,
            self.start,
            self.ref,
            self.alt,
            self.reference_name)

    def __repr__(self):
        return str(self)

    def __hash__(self):
        return self.start

    def __lt__(self, other: object) -> bool:
        """
        Variants are ordered by locus.
        """
        if not isinstance(other, Variant):
            return NotImplemented
        if self.contig == other.contig:
            return self.start < other.start
        return self.contig < other.contig

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Variant):
            return NotImplemented
        if self is other:
            return True
        return (
            self.contig == other.contig and
            self.start == other.start and
            self.end == other.end and
            self.ref == other.ref and
            self.alt == other.alt and
            self.genome == other.genome)

    def to_dict(self):
        """
        We want the original values (un-normalized) field values while
        serializing since normalization will happen in __init__.
        """
        state = dict(
            contig=self.original_contig,
            start=self.original_start,
            ref=self.original_ref,
            alt=self.original_alt,
            genome=self.genome,
            original_reference_name=self.original_reference_name,
            original_genome_was_ucsc=self.original_genome_was_ucsc,
            original_genome=(self.original_genome if isinstance(self.original_genome, (str, int)) else None),
            allow_extended_nucleotides=self.allow_extended_nucleotides,
            normalize_contig_names=self.normalize_contig_names,
            convert_ucsc_contig_names=self.convert_ucsc_contig_names)
        if self.chrm_is_rcrs:
            state["chrm_is_rcrs"] = True
        return state

    @property
    def trimmed_ref(self):
        """
        Eventually the field Variant.ref will store the reference nucleotides
        as given in a VCF or MAF and trimming of any shared prefix/suffix
        between ref and alt will be done via the properties trimmed_ref
        and trimmed_alt.
        """
        return self.ref

    @property
    def trimmed_alt(self):
        """
        Eventually the field Variant.ref will store the reference nucleotides
        as given in a VCF or MAF and trimming of any shared prefix/suffix
        between ref and alt will be done via the properties trimmed_ref
        and trimmed_alt.
        """
        return self.alt

    @property
    def trimmed_base1_start(self):
        """
        Currently the field Variant.start carries the base-1 starting position
        adjusted by trimming any shared prefix between Variant.ref and
        Variant.alt. Eventually this trimming should be done more explicitly
        via trimmed_* properties.
        """
        return self.start

    @property
    def trimmed_base1_end(self):
        """
        Currently the field Variant.end carries the base-1 "last" position of
        this variant, adjusted by trimming any shared suffix between
        Variant.ref and Variant.alt. Eventually this trimming should be done
        more explicitly via trimmed_* properties.
        """
        return self.end

    @property
    def short_description(self):
        """
        HGVS nomenclature for genomic variants
        More info: http://www.hgvs.org/mutnomen/
        """
        if self.is_insertion:
            return "chr%s g.%d_%dins%s" % (
                self.contig,
                self.start,
                self.start + 1,
                self.alt)
        elif self.is_deletion:
            # HGVS uses a single coordinate for a one-base deletion.
            # https://hgvs-nomenclature.org/stable/recommendations/DNA/deletion/
            if self.start == self.end:
                return "chr%s g.%ddel%s" % (self.contig, self.start, self.ref)
            return "chr%s g.%d_%ddel%s" % (
                self.contig,
                self.start,
                self.end,
                self.ref)
        elif self.ref == self.alt:
            # Equal alleles trim to empty strings. Describe the original
            # checked interval, using HGVS's explicit no-change marker.
            # https://hgvs-nomenclature.org/stable/recommendations/DNA/substitution/
            end = self.original_start + max(0, len(self.original_ref) - 1)
            interval = str(self.original_start)
            if end != self.original_start:
                interval += "_%d" % end
            return "chr%s g.%s=" % (self.contig, interval)
        else:
            # substitution
            return "chr%s g.%d%s>%s" % (
                self.contig,
                self.start,
                self.ref,
                self.alt)


    def _check_that_genome_has_contig(self):
        """
        Any annotation code which wants to look up genomic loci in a reference
        genome should first run this helper to make sure the normalized contig
        name is in the genome.
        """
        # Different annotation datasets can share an assembly name. Query
        # this genome's contigs; pyensembl caches the query per database.
        if self.contig not in self.genome.contigs():
            message = "Invalid contig name '%s' for reference '%s'" % (
                self.contig, self.reference_name)
            note = self._missing_contig_note()
            if note:
                message += ": " + note
            raise InvalidContigError(message)

    def _missing_contig_note(self):
        """Why this variant's mitochondrion has no annotation, or None.

        hg19's chrM is not the rCRS (#707). T2T-CHM13v2.0's own mitochondrion
        is unannotated; only the rCRS of its ``_rCRS`` analysis set is (#722).
        """
        assembly = getattr(self.genome, "assembly", self.reference_name)
        if _is_hg19_mitochondrion(self.contig, assembly, self.chrm_is_rcrs):
            return _HG19_MITOCHONDRION_NOTE
        return mitochondrion_note(self.contig, assembly)



    @property
    def transcripts(self):
        """
        Return Transcript objects for all transcripts which overlap this
        variant. A transcript joined across a circular origin does not
        overlap the stretch between its 3' and 5' ends.
        """
        if self._transcripts is None:
            self._check_that_genome_has_contig()
            self._transcripts = circular_origin.transcripts_at_locus(
                self.genome, self.contig, self.start, self.end)
        return self._transcripts

    @property
    def coding_transcripts(self):
        """
        Protein coding transcripts
        """
        return [
            transcript
            for transcript in self.transcripts
            if transcript.is_protein_coding
        ]

    @property
    def transcript_ids(self):
        return [transcript.id for transcript in self.transcripts]

    @property
    def transcript_names(self):
        return [transcript.name for transcript in self.transcripts]

    @property
    def genes(self):
        """
        Return Gene object for all genes which overlap this variant.
        """
        if self._genes is None:
            self._check_that_genome_has_contig()
            self._genes = circular_origin.genes_at_locus(
                self.genome, self.contig, self.start, self.end)
        return self._genes

    @property
    def gene_ids(self):
        """
        Return IDs of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        if self._gene_ids is None:
            self._check_that_genome_has_contig()
            self._gene_ids = circular_origin.gene_ids_at_locus(
                self.genome, self.contig, self.start, self.end)
        return self._gene_ids

    @property
    def gene_names(self):
        """
        Return names of all genes which overlap this variant. Calling
        this method is significantly cheaper than calling `Variant.genes()`,
        which has to issue many more queries to construct each Gene object.
        """
        if self._gene_names is None:
            self._check_that_genome_has_contig()
            self._gene_names = circular_origin.gene_names_at_locus(
                self.genome, self.contig, self.start, self.end)
        return self._gene_names

    @property
    def coding_genes(self):
        """
        Protein coding transcripts
        """
        return [
            gene for gene in self.genes
            if gene.is_protein_coding
        ]

    def orf_effects(self, *, catalog=None, transcripts=None, **kwargs):
        """Additive reference-ORF consequences; see :mod:`varcode.orf`.

        Reference MS observations never imply mutant translation.
        """
        from .orf import predict_orf_effects
        return predict_orf_effects(self, catalog=catalog, transcripts=transcripts, **kwargs)

    def effects(
            self, raise_on_error: bool = True,
            annotator: str | EffectAnnotator | None = None,
            phase_resolver=None, rna_resolver=None,
            germline=None, extra_orfs=None, transcript_filter=None,
    ) -> EffectCollection:
        """Predict the variant's effects on overlapping transcripts.

        The default composes phase, splice and local structural edits into
        baseline/mutant transcript comparisons. The ranked effect retains
        alternatives in ``candidates``. Explicit ``annotator="fast"`` preserves
        the legacy splice/phase/SV wrappers.

        Parameters
        ----------
        raise_on_error : bool
            If True, raise on annotation errors; if False, capture
            per-transcript errors as Failure effects. Failed initial
            gene/transcript lookups return a transcript-free Failure
            and log the error.

        annotator : str, EffectAnnotator, or None
            Per-call annotator override. ``None`` uses the currently
            configured default (``"transcript_model"``; swappable via
            :func:`varcode.set_default_annotator` or
            :func:`varcode.use_annotator`). String names are resolved
            against the registry. See openvax/varcode#271.

        phase_resolver : resolver object or None
            Optional source of cis/trans constraints and externally supplied
            transcript models. Observed models refine the same candidates;
            partial/unmatched models stay unresolved with their actual sequence.
            RNA sequence does not establish translation or missing phase.

        rna_resolver : RNAEvidenceResolver or None
            Optional source of observed outcomes. Exact full-cDNA matches
            support compatible hypotheses without discarding other baselines.
            Missing sequence remains unavailable; counts are not probabilities.
            Legacy annotators retain their documented reconciliation behavior.

        germline : GermlineContext or None
            Supplied patient alleles across the transcript span. Unknown phase
            produces candidate comparisons or HypothesisLimit; known cis alleles
            modify both baseline and mutant. An identical inherited allele
            returns GermlineAlleleOverlap, with LOH unassessed. Omission uses the
            annotation reference and cannot recover missing patient information.

        extra_orfs : str, ORFdata rule, ORFCatalog or iterable of ReferenceORF, optional
            Additional reference coding frames, supplied by ORFdata. ``"all"``
            uses every ORF ORFdata has for this variant's genome: its species
            and assembly's catalogues, lifted from another assembly when
            needed, placed on every transcript carrying each ORF's exact coding
            path. An evidence name selects among them: ``"ms"`` (any MS
            peptide), ``"ribosome"`` or ``"riboseq"`` (Ribo-seq),
            ``"proteomics"`` (peptides from digested cellular proteins),
            ``"immunopeptidomics"`` (peptides eluted from MHC molecules) or
            ``"ms_unclassified"`` (MS whose source doesn't say which). A list of names, such as
            ``["proteomics", "immunopeptidomics"]``, selects ORFs with any of
            them, and a rule such as ``orfdata.RiboSeq() & orfdata.MS()``
            those meeting all of its parts. Placements are made once per genome and reused. A catalogue
            is resolved against this variant's genome; bound ORFs retain their
            transcript placements. Ordinary effects are preserved, and native
            ORF consequences are appended using the transcript-model engine.
            Each has ``orf``, ``reference_transcript`` and ``orf_annotation``
            attributes. No catalogue is loaded when omitted. ``annotator``
            selects the ordinary transcript annotator; additional ORFs use the
            existing ORF engine, including germline and DNA phase hypotheses.
            Observed RNA requires ORF-specific base mapping through
            ``annotate_orf_rna_product``; incompatible RNA observations raise
            or become Failure effects according to ``raise_on_error``.

        transcript_filter : callable or None
            Select transcripts before annotation, for example
            ``lambda t: t.biotype == "protein_coding"`` or
            ``lambda t: t.id in selected_ids``. Omission keeps all overlapping
            transcripts. The predicate receives the original reference
            transcript, including for additional ORF placements. Selection
            does not change coding-biotype or completeness rules. Genuine
            transcript-free overlaps and unplaced ORF failures remain visible;
            rejecting all overlapping transcripts can return an empty collection.
            Predicate errors follow ``raise_on_error``. Use a deterministic
            predicate: it may be called again for another variant or ORF.

        Returns
        -------
        EffectCollection
            Transcript consequences and any requested additional ORF effects.
            A consequence may retain several candidate outcomes; ranking does
            not establish which one is expressed. Transcript selection can
            produce an empty collection.
        """
        from .extra_orfs import prepare_extra_orfs
        prepared = prepare_extra_orfs(extra_orfs)  # A bad value raises before any annotation.
        effects = predict_variant_effects(
            variant=self,
            raise_on_error=raise_on_error,
            annotator=annotator,
            germline=germline,
            phase_resolver=phase_resolver,
            transcript_filter=transcript_filter,
        )
        if phase_resolver is not None:
            from .phasing import apply_phase_resolver_to_effects
            apply_phase_resolver_to_effects(effects, phase_resolver)
        if rna_resolver is not None:
            from .rna_evidence import apply_rna_evidence_to_effects
            apply_rna_evidence_to_effects(effects, rna_resolver)
        if prepared is not None:
            additional = prepared.effects(
                self, raise_on_error=raise_on_error, germline=germline,
                phase_resolver=phase_resolver, rna_resolver=rna_resolver,
                transcript_filter=transcript_filter)
            effects = effects.clone_with_new_elements([*effects, *additional])
        return effects

    def effect_on_transcript(
            self, transcript: Transcript,
            annotator: str | EffectAnnotator | None = None,
            germline=None, phase_resolver=None,
    ) -> MutationEffect:
        """Annotate one transcript using the same selection as :meth:`effects`.

        Parameters
        ----------
        transcript : pyensembl.Transcript
            Reference transcript to annotate.
        annotator : str, EffectAnnotator, or None
            Registered name or instance; None uses the current default.
        germline : GermlineContext or None
            Optional patient context for the selected annotator.
        phase_resolver : resolver object or None
            Optional evidence for germline phase.

        Returns
        -------
        MutationEffect
            Consequence, including an outcome set where needed. An annotator's
            explicit refusal becomes Unresolved; it does not return None or
            NotImplemented from this public method.

        Notes
        -----
        Annotation errors propagate. Use ``effects(raise_on_error=False)``
        to collect Failure effects across overlapping transcripts instead.
        """
        return predict_variant_effect_on_transcript(
            self, transcript, annotator=annotator, germline=germline,
            phase_resolver=phase_resolver)

    @property
    def is_insertion(self):
        """
        Does this variant represent the insertion of nucleotides into the
        reference genome?
        """
        # An insertion would appear in a VCF like C>CT, so that the
        # alternate allele starts with the reference nucleotides.
        # Since the nucleotide strings may be normalized in the constructor,
        # it's worth noting that the normalized form of this variant would be
        # ''>'T', so that 'T'.startswith('') still holds.
        return (len(self.ref) < len(self.alt)) and self.alt.startswith(self.ref)

    @property
    def is_deletion(self):
        """
        Does this variant represent the deletion of nucleotides from the
        reference genome?
        """
        # A deletion would appear in a VCF like CT>C, so that the
        # reference allele starts with the alternate nucleotides.
        # This is true even in the normalized case, where the alternate
        # nucleotides are an empty string.
        return (len(self.ref) > len(self.alt)) and self.ref.startswith(self.alt)

    @property
    def is_indel(self):
        """Is this variant either an insertion or deletion?"""
        return self.is_insertion or self.is_deletion

    @property
    def is_snv(self):
        """Is the variant a single nucleotide variant"""
        return (len(self.ref) == len(self.alt) == 1) and (self.ref != self.alt)

    @property
    def is_transition(self):
        """Is this variant and pyrimidine to pyrimidine change or purine to purine change"""
        return self.is_snv and is_purine(self.ref) == is_purine(self.alt)

    @property
    def is_transversion(self):
        """Is this variant a pyrimidine to purine change or vice versa"""
        return self.is_snv and is_purine(self.ref) != is_purine(self.alt)

    def _normalize_contig_name(self, contig):
        """
        Apply this variant's contig-name settings to ``contig``: tidy the
        name when ``normalize_contig_names`` is set, then trim the "chr"
        prefix from UCSC names to match Ensembl (and convert "chrM" to
        "MT") when ``convert_ucsc_contig_names`` is set. Shared with
        StructuralVariant so breakend mates are named like the variant's
        own contig.

        Under GRCh37, conversion leaves ``chrM`` and ``M`` as ``chrM``:
        hg19's chrM is not the rCRS that Ensembl names MT, and its positions
        are offset from it (#707). Such a variant fails annotation, unless
        ``chrm_is_rcrs`` declares this chrM the rCRS (#710).
        """
        assembly = getattr(self.genome, "assembly", self.reference_name)
        if self.normalize_contig_names:
            return _annotation_contig_name(
                contig, assembly, self.convert_ucsc_contig_names, self.chrm_is_rcrs)
        if not self.convert_ucsc_contig_names:
            return contig
        pinned = pinned_contig_name(contig, assembly, self.chrm_is_rcrs)
        if pinned is not None:
            return pinned
        if _is_hg19_mitochondrion(contig, assembly, self.chrm_is_rcrs):
            return "chrM"
        return self._convert_ucsc_contig_name_to_ensembl(contig)

    def _convert_ucsc_contig_name_to_ensembl(self, contig):
        """
        Convert UCSC chromosome and accession-bearing scaffold names.

        Preserve explicit accession versions. The unversioned GL names in
        hg19 are version 1 (UCSC hg19.chromAlias.txt); do not reinterpret them
        as a different version found in a newer assembly. This is a name
        conversion only and does not open or download reference data.

        Parameters
        ----------
        contig : str

        Returns
        -------
        str
        """
        return _convert_ucsc_contig_name(
            contig, getattr(self.genome, "assembly", self.reference_name))

    def clone_without_ucsc_data(self):
        """
        Clone this variant but discarding the original format of its genome
        and contig: useful if we want to mix hg19 and GRCh37 variants.

        Returns
        -------
        Variant
        """
        return Variant(
            contig=self.contig,
            start=self.original_start,
            ref=self.original_ref,
            alt=self.original_alt,
            genome=self.genome,
            allow_extended_nucleotides=self.allow_extended_nucleotides,
            normalize_contig_names=self.normalize_contig_names,
            convert_ucsc_contig_names=False)

def variant_ascending_position_sort_key(variant):
    """
    Sort key function used to sort variants in ascending order by
    chromosomal position.

    This is a function rather than method to make it easier to serialize.

    Returns
    -------
    tuple of (str, int)
    """
    return (variant.contig, variant.start)
