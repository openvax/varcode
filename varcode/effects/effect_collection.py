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

from collections import OrderedDict

import pandas as pd
from sercol import Collection

from ..csv_helpers import (
    CONTIG_COLUMN_ALIASES,
    read_metadata_header,
    resolve_contig_column,
    warn_on_version_drift,
    write_metadata_header,
)
from ..version import __version__ as _varcode_version
from .effect_ordering import (
    effect_priority,
    multi_gene_effect_sort_key,
    top_priority_effect,
    transcript_effect_priority_dict
)


def _default_effect_sort_key(effect):
    """Default sort key for EffectCollection: highest-priority effects
    first. Negating the priority integer turns Python's ascending sort
    into descending-by-priority. See openvax/varcode#227.
    """
    return -effect_priority(effect)


class EffectCollection(Collection):
    """
    Collection of MutationEffect objects and helpers for grouping or filtering
    them.
    """
    def __init__(
            self,
            effects,
            distinct=False,
            sort_key=None,
            sources=set([]),
            annotator=None,
            annotator_version=None,
            annotated_at=None):
        """
        Parameters
        ----------
        effects : list
            Collection of any class which  is compatible with the sort key


        distinct : bool
            Only keep distinct entries or allow duplicates.

        sort_key : fn
            Function which maps each element to a sorting criterion.
            If None (the default), effects are sorted by priority with
            the most severe effects first. Pass an explicit sort_key to
            override this behaviour, or `False` to disable sorting.

        sources : set
            Set of files from which this collection was generated.

        annotator : str or None
            Name of the :class:`EffectAnnotator` that produced these
            effects. Populated automatically by
            :func:`predict_variant_effects`; ``None`` for collections
            built by hand. See openvax/varcode#271.

        annotator_version : str or None
            Version string of the annotator (typically the varcode
            version for built-in annotators). ``None`` when
            ``annotator`` is None.

        annotated_at : str or None
            ISO-8601 UTC timestamp recording when the annotation ran.
            Populated by :func:`predict_variant_effects`.
        """
        if sort_key is None:
            sort_key = _default_effect_sort_key
        elif sort_key is False:
            sort_key = None
        Collection.__init__(
            self,
            elements=effects,
            distinct=distinct,
            sort_key=sort_key,
            sources=sources)
        # Keep self.effects in sync with the Collection's post-sort
        # elements so that iterating and reading `.effects` produce
        # the same order.  See openvax/varcode#220.
        self.effects = self.elements
        self.annotator = annotator
        self.annotator_version = annotator_version
        self.annotated_at = annotated_at

    def to_dict(self):
        return dict(
            effects=self.effects,
            sort_key=self.sort_key,
            distinct=self.distinct,
            sources=self.sources,
            annotator=self.annotator,
            annotator_version=self.annotator_version,
            annotated_at=self.annotated_at)

    def clone_with_new_elements(self, new_elements):
        # Forward annotator provenance explicitly — sercol's base
        # clone_with_new_elements only carries distinct/sort_key/sources,
        # so without this the cloned collection would silently lose
        # its #271 provenance metadata.
        return Collection.clone_with_new_elements(
            self,
            new_elements,
            rename_dict={"elements": "effects"},
            extra_kwargs={
                "annotator": self.annotator,
                "annotator_version": self.annotator_version,
                "annotated_at": self.annotated_at,
            })

    def groupby_variant(self):
        return self.groupby(key_fn=lambda effect: effect.variant)

    def groupby_transcript(self):
        return self.groupby(key_fn=lambda effect: effect.transcript)

    def groupby_transcript_name(self):
        return self.groupby(key_fn=lambda effect: effect.transcript_name)

    def groupby_transcript_id(self):
        return self.groupby(key_fn=lambda effect: effect.transcript_id)

    def groupby_gene(self):
        return self.groupby(key_fn=lambda effect: effect.gene)

    def groupby_gene_name(self):
        return self.groupby(key_fn=lambda effect: effect.gene_name)

    def groupby_gene_id(self):
        return self.groupby(key_fn=lambda effect: effect.gene_id)

    def gene_counts(self):
        """
        Returns number of elements overlapping each gene name. Expects the
        derived class (VariantCollection or EffectCollection) to have
        an implementation of groupby_gene_name.
        """
        return {
            gene_name: len(group)
            for (gene_name, group)
            in self.groupby_gene_name().items()
        }

    def filter_by_transcript_expression(
            self,
            transcript_expression_dict,
            min_expression_value=0.0):
        """
        Filters effects to those which have an associated transcript whose
        expression value in the transcript_expression_dict argument is greater
        than min_expression_value.

        Parameters
        ----------
        transcript_expression_dict : dict
            Dictionary mapping Ensembl transcript IDs to expression estimates
            (either FPKM or TPM)

        min_expression_value : float
            Threshold above which we'll keep an effect in the result collection
        """
        return self.filter_above_threshold(
            key_fn=lambda effect: effect.transcript_id,
            value_dict=transcript_expression_dict,
            threshold=min_expression_value)

    def filter_by_gene_expression(
            self,
            gene_expression_dict,
            min_expression_value=0.0):
        """
        Filters effects to those which have an associated gene whose
        expression value in the gene_expression_dict argument is greater
        than min_expression_value.

        Parameters
        ----------
        gene_expression_dict : dict
            Dictionary mapping Ensembl gene IDs to expression estimates
            (either FPKM or TPM)

        min_expression_value : float
            Threshold above which we'll keep an effect in the result collection
        """
        return self.filter_above_threshold(
            key_fn=lambda effect: effect.gene_id,
            value_dict=gene_expression_dict,
            threshold=min_expression_value)

    def filter_by_effect_priority(self, min_priority_class):
        """
        Create a new EffectCollection containing only effects whose priority
        falls below the given class.
        """
        min_priority = transcript_effect_priority_dict[min_priority_class]
        return self.filter(
            lambda effect: effect_priority(effect) >= min_priority)

    def drop_silent_and_noncoding(self):
        """
        Create a new EffectCollection containing only non-silent coding effects
        """
        return self.filter(lambda effect: effect.modifies_protein_sequence)

    def detailed_string(self):
        """
        Create a long string with all transcript effects for each mutation,
        grouped by gene (if a mutation affects multiple genes).
        """
        lines = []
        # TODO: annoying to always write `groupby_result.items()`,
        # consider makings a GroupBy class which iterates over pairs
        # and also common helper methods like `map_values`.
        for variant, variant_effects in self.groupby_variant().items():
            lines.append("\n%s" % variant)

            gene_effects_groups = variant_effects.groupby_gene_id()
            for (gene_id, gene_effects) in gene_effects_groups.items():
                if gene_id:
                    gene_name = variant.ensembl.gene_name_of_gene_id(gene_id)
                    lines.append("  Gene: %s (%s)" % (gene_name, gene_id))
                # place transcript effects with more significant impact
                # on top (e.g. FrameShift should go before NoncodingTranscript)
                for effect in sorted(
                        gene_effects,
                        key=effect_priority,
                        reverse=True):
                    lines.append("  -- %s" % effect)

            # if we only printed one effect for this gene then
            # it's redundant to print it again as the highest priority effect
            if len(variant_effects) > 1:
                best = variant_effects.top_priority_effect()
                lines.append("  Highest Priority Effect: %s" % best)
        return "\n".join(lines)

    def top_priority_effect(self):
        """Highest priority MutationEffect of all genes/transcripts overlapped
        by this variant. If this variant doesn't overlap anything, then this
        this method will return an Intergenic effect.

        If multiple effects have the same priority, then return the one
        which is associated with the longest transcript.
        """
        return top_priority_effect(self.elements)

    # TODO: find a way to express these kinds of methods without
    # duplicating every single groupby_* method

    def top_priority_effect_per_variant(self):
        """Highest priority effect for each unique variant"""
        return OrderedDict(
            (variant, top_priority_effect(variant_effects))
            for (variant, variant_effects)
            in self.groupby_variant().items())

    def top_priority_effect_per_transcript_id(self):
        """Highest priority effect for each unique transcript ID"""
        return OrderedDict(
            (transcript_id, top_priority_effect(variant_effects))
            for (transcript_id, variant_effects)
            in self.groupby_transcript_id().items())

    def top_priority_effect_per_gene_id(self):
        """Highest priority effect for each unique gene ID"""
        return OrderedDict(
            (gene_id, top_priority_effect(variant_effects))
            for (gene_id, variant_effects)
            in self.groupby_gene_id().items())

    def effect_expression(self, expression_levels):
        """
        Parameters
        ----------
        expression_levels : dict
            Dictionary mapping transcript IDs to length-normalized expression
            levels (either FPKM or TPM)

        Returns
        -------
        OrderedDict
            Mapping from each transcript effect to an expression quantity.
            Effects that don't have an associated transcript (e.g. Intergenic)
            are excluded.
        """
        return OrderedDict(
            (effect, expression_levels.get(effect.transcript.id, 0.0))
            for effect in self
            if effect.transcript is not None)

    def top_expression_effect(self, expression_levels):
        """
        Return effect whose transcript has the highest expression level.
        If none of the effects are expressed or have associated transcripts,
        then return None. In case of ties, add lexicographical sorting by
        effect priority and transcript length.
        """
        effect_expression_dict = self.effect_expression(expression_levels)

        if len(effect_expression_dict) == 0:
            return None

        def key_fn(effect_fpkm_pair):
            """
            Sort effects primarily by their expression level
            and secondarily by the priority logic used in
            `top_priority_effect`.
            """
            (effect, fpkm) = effect_fpkm_pair
            return (fpkm, multi_gene_effect_sort_key(effect))

        return max(effect_expression_dict.items(), key=key_fn)[0]

    _DATAFRAME_COLUMNS = (
        "variant",
        "contig", "start", "ref", "alt",
        "is_snv", "is_transversion", "is_transition",
        "gene_id", "gene_name",
        "transcript_id", "transcript_name",
        "effect_type", "effect",
    )

    def to_dataframe(self):
        """Build a dataframe from the effect collection."""
        # list of properties to extract from Variant objects if they're
        # not None
        variant_properties = [
            "contig",
            "start",
            "ref",
            "alt",
            "is_snv",
            "is_transversion",
            "is_transition"
        ]

        def row_from_effect(effect):
            row = OrderedDict()

            row['variant'] = str(effect.variant.short_description)
            for field_name in variant_properties:
                # if effect.variant is None then this column value will be None
                row[field_name] = getattr(effect.variant, field_name, None)
            row['gene_id'] = effect.gene_id
            row['gene_name'] = effect.gene_name
            row['transcript_id'] = effect.transcript_id
            row['transcript_name'] = effect.transcript_name

            row['effect_type'] = effect.__class__.__name__
            row['effect'] = effect.short_description
            return row
        # Always emit the same column set even for empty collections so
        # CSV round-trip doesn't have to special-case len == 0.
        return pd.DataFrame.from_records(
            [row_from_effect(effect) for effect in self],
            columns=self._DATAFRAME_COLUMNS,
        )

    def to_csv(self, path, include_header=True):
        """Write this collection to CSV.

        Parameters
        ----------
        path : str
            Output path.

        include_header : bool
            If True (default), prepend ``# key=value`` metadata lines
            with varcode version and reference genome so the file can
            be read back via ``from_csv`` without supplying a genome
            explicitly. Pass ``include_header=False`` for fast
            consumers that don't tolerate comment lines.
        """
        df = self.to_dataframe()
        if not include_header:
            df.to_csv(path, index=False)
            return

        metadata = OrderedDict()
        metadata["varcode_version"] = _varcode_version
        metadata["reference_name"] = self._serialized_reference_name()
        # Annotator provenance (#271 stage 3b). Each field is skipped
        # when None by write_metadata_header.
        metadata["annotator"] = self.annotator
        metadata["annotator_version"] = self.annotator_version
        metadata["annotated_at"] = self.annotated_at
        with open(path, "w") as f:
            write_metadata_header(f, metadata)
            df.to_csv(f, index=False)

    def _serialized_reference_name(self):
        """Pick a reference-name string to record in the CSV header.

        Uses the first effect's variant's reference_name. Returns None
        if the collection is empty or has no consistent reference.
        """
        names = set()
        for effect in self:
            variant = getattr(effect, "variant", None)
            if variant is not None and variant.reference_name:
                names.add(variant.reference_name)
        if len(names) == 1:
            return next(iter(names))
        return None

    @classmethod
    def from_csv(cls, path, genome=None):
        """Rebuild an EffectCollection from a CSV previously written by
        ``EffectCollection.to_csv()``.

        The current CSV format records (contig, start, ref, alt,
        transcript_id) but not enough per-effect state to reconstruct
        effects byte-for-byte. This method takes the pragmatic semantic
        round-trip path: rebuild each Variant, re-annotate against the
        recorded transcript, and emit the resulting effect. The
        resulting collection should match the original whenever
        annotation is deterministic for a given (variant, transcript)
        pair.

        **Prefer ``from_json`` for byte-for-byte round-trip** or for
        larger collections (≳10k effects); per-row re-annotation makes
        CSV loading significantly slower than JSON. Emits a warning
        when the CSV header reports a different *major* varcode version
        than the one currently installed — annotation logic can change
        across major versions and the reconstructed effects may differ
        from the ones that were written.

        Parameters
        ----------
        path : str
            Path to the CSV file. Lines starting with '#' are treated as
            comments and parsed as ``key=value`` metadata.

        genome : pyensembl.Genome, str, int, or None
            Reference genome to associate with the loaded variants and
            to look up transcripts by ID. If ``None``, the reference is
            read from the CSV's metadata header
            (``# reference_name=...``). If neither is available, raises
            ``ValueError``.

        Returns
        -------
        EffectCollection
        """
        # Import here to avoid a circular import at module load time.
        import warnings
        from ..variant import Variant

        header = read_metadata_header(path)
        warn_on_version_drift(header, _varcode_version, path)
        if genome is None:
            genome = header.get("reference_name")
        if genome is None:
            raise ValueError(
                "from_csv needs a reference genome: pass the `genome` "
                "argument explicitly, or write the CSV with "
                "`to_csv(include_header=True)` so `# reference_name=...` "
                "is recorded in the header. Neither was found at %s." % path)

        # Accept either "contig" or "chr" as the contig column so
        # CSVs are interchangeable between VariantCollection and
        # EffectCollection (openvax/varcode#274). Declaring both in
        # dtype is a no-op for whichever column is absent.
        df = pd.read_csv(
            path, comment="#", dtype={"chr": str, "contig": str})
        contig_col = resolve_contig_column(df.columns)
        if contig_col is None:
            raise ValueError(
                "CSV at %s is missing a contig column: expected one of %s."
                % (path, list(CONTIG_COLUMN_ALIASES)))
        required = {"start", "ref", "alt", "transcript_id", "effect_type"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                "CSV at %s is missing required columns: %s" % (
                    path, sorted(missing)))

        # Extract columns by name and coerce types up front; zip is
        # robust against column reordering and extra columns, unlike
        # itertuples which depends on attribute access to
        # valid-identifier column names in fixed positions.
        contigs = df[contig_col].astype(str)
        starts = df["start"].astype(int)
        refs = df["ref"].fillna("")
        alts = df["alt"].fillna("")
        transcript_ids = df["transcript_id"]
        effect_types = df["effect_type"]

        effects = []
        resolved_genome = None
        for contig, start, ref, alt, transcript_id, effect_type in zip(
                contigs, starts, refs, alts, transcript_ids, effect_types):
            variant = Variant(
                contig=contig,
                start=start,
                ref=ref,
                alt=alt,
                genome=genome,
            )
            # Cache the resolved pyensembl.Genome from the first variant
            # so subsequent transcript lookups don't pay the inference cost.
            if resolved_genome is None:
                resolved_genome = variant.ensembl
            if pd.isna(transcript_id) or transcript_id == "":
                # Row is for an intergenic / intragenic effect — no
                # transcript context. Rebuild by running the full effect
                # set on the variant and taking the non-transcript match.
                effects_for_variant = variant.effects()
                matched = False
                for e in effects_for_variant:
                    if e.transcript is None and e.__class__.__name__ == effect_type:
                        effects.append(e)
                        matched = True
                        break
                if not matched:
                    warnings.warn(
                        "Could not recover effect of type %r for variant %s "
                        "with no transcript_id in %s; row dropped from "
                        "reconstructed collection." % (
                            effect_type, variant, path)
                    )
                continue
            transcript = resolved_genome.transcript_by_id(str(transcript_id))
            effects.append(variant.effect_on_transcript(transcript))
        return cls(
            effects=effects,
            annotator=header.get("annotator"),
            annotator_version=header.get("annotator_version"),
            annotated_at=header.get("annotated_at"),
        )
