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

from .effects import EffectCollection
from .common import memoize
from .csv_helpers import (
    CONTIG_COLUMN_ALIASES,
    read_metadata_header,
    resolve_contig_column,
    write_metadata_header,
)
from .variant import Variant, variant_ascending_position_sort_key
from .version import __version__ as _varcode_version


class VariantCollection(Collection):
    def __init__(
            self,
            variants,
            distinct=True,
            sort_key=variant_ascending_position_sort_key,
            sources=None,
            source_to_metadata_dict={}):
        """
        Construct a VariantCollection from a list of Variant records.

        Parameters
        ----------
        variants : iterable
            Variant objects contained in this VariantCollection

        distinct : bool
            Don't keep repeated variants

        sort_key : callable

        sources : set
            Optional set of source names, may be larger than those for
            which we have metadata dictionaries.

        source_to_metadata_dict : dict
            Dictionary mapping each source name (e.g. VCF path) to a dictionary
            from metadata attributes to values.
        """
        self.source_to_metadata_dict = source_to_metadata_dict
        if sources is None:
            sources = set(source_to_metadata_dict.keys())
        if any(source not in sources for source in source_to_metadata_dict.keys()):
            raise ValueError(
                "Mismatch between sources=%s and keys of source_to_metadata_dict=%s" % (
                    sources,
                    set(source_to_metadata_dict.keys())))
        Collection.__init__(
            self,
            elements=variants,
            distinct=distinct,
            sort_key=sort_key,
            sources=sources)
        # Keep self.variants in sync with the Collection's post-sort,
        # post-dedup elements so that iterating `vc` and reading
        # `vc.variants` produce the same order.  See openvax/varcode#220.
        self.variants = self.elements

    @property
    def metadata(self):
        """
        The most common usage of a VariantCollection is loading a single VCF,
        in which case it's annoying to have to always specify that path
        when accessing metadata fields. This property is meant to both maintain
        backward compatibility with old versions of Varcode and make the common
        case easier.
        """
        return self.source_to_metadata_dict[self.source]

    def to_dict(self):
        """
        Since Collection.to_dict() returns a state dictionary with an
        'elements' field we have to rename it to 'variants'.
        """
        return dict(
            variants=self.variants,
            distinct=self.distinct,
            sort_key=self.sort_key,
            sources=self.sources,
            source_to_metadata_dict=self.source_to_metadata_dict)

    def clone_with_new_elements(self, new_elements):
        """
        Create another VariantCollection of the same class and with
        same state (including metadata) but possibly different entries.

        Warning: metadata is a dictionary keyed by variants. This method
        leaves that dictionary as-is, which may result in extraneous entries
        or missing entries.
        """
        kwargs = self.to_dict()
        kwargs["variants"] = new_elements
        return self.from_dict(kwargs)

    def effects(self, raise_on_error=True):
        """
        Parameters
        ----------
        raise_on_error : bool, optional
            If exception is raised while determining effect of variant on a
            transcript, should it be raised? This default is True, meaning
            errors result in raised exceptions, otherwise they are only logged.

        """
        return EffectCollection([
            effect
            for variant in self
            for effect in variant.effects(raise_on_error=raise_on_error)
        ])

    @memoize
    def reference_names(self):
        """
        All distinct reference names used by Variants in this
        collection.

        Returns
        -------
        set of str
        """
        return {variant.reference_name for variant in self}

    @memoize
    def original_reference_names(self):
        """
        Similar to `reference_names` but preserves UCSC references,
        so that a variant collection derived from an hg19 VCF would
        return {"hg19"} instead of {"GRCh37"}.

        Returns
        -------
        set of str
        """
        return {variant.original_reference_name for variant in self}

    def groupby_gene(self):
        return self.multi_groupby(lambda x: x.genes)

    def groupby_gene_name(self):
        """
        Group variants by the gene names they overlap, which may put each
        variant in multiple groups.
        """
        return self.multi_groupby(lambda x: x.gene_names)

    def groupby_gene_id(self):
        return self.multi_groupby(lambda x: x.gene_ids)

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

    def detailed_string(self):
        lines = []
        gene_groups = self.groupby_gene_id()
        for gene_id, variant_group in sorted(gene_groups.items()):
            # in case we are combining variants with different Ensembl releases
            # use them all to gather gene names
            gene_names = {
                variant.ensembl.gene_name_of_gene_id(gene_id)
                for variant in variant_group
            }
            gene_name_str = ", ".join(sorted(gene_names))
            lines.append("  -- %s (%s):" % (gene_name_str, gene_id))
            for variant in variant_group:
                lines.append("  ---- %s" % variant)
        header = self.short_string()
        joined_lines = "\n".join(lines)
        return "%s\n%s" % (header, joined_lines)

    def filter_by_transcript_expression(
            self,
            transcript_expression_dict,
            min_expression_value=0.0):
        """
        Filters variants down to those which have overlap a transcript whose
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
        return self.filter_any_above_threshold(
            multi_key_fn=lambda variant: variant.transcript_ids,
            value_dict=transcript_expression_dict,
            threshold=min_expression_value)

    def filter_by_gene_expression(
            self,
            gene_expression_dict,
            min_expression_value=0.0):
        """
        Filters variants down to those which have overlap a gene whose
        expression value in the transcript_expression_dict argument is greater
        than min_expression_value.

        Parameters
        ----------
        gene_expression_dict : dict
            Dictionary mapping Ensembl gene IDs to expression estimates
            (either FPKM or TPM)

        min_expression_value : float
            Threshold above which we'll keep an effect in the result collection
        """
        return self.filter_any_above_threshold(
            multi_key_fn=lambda effect: effect.gene_ids,
            value_dict=gene_expression_dict,
            threshold=min_expression_value)

    def exactly_equal(self, other):
        '''
        Comparison between VariantCollection instances that takes into account
        the info field of Variant instances.

        Returns
        ----------
        True if the variants in this collection equal the variants in the other
        collection. The Variant.info fields are included in the comparison.
        '''
        return (
            self.__class__ == other.__class__ and
            len(self) == len(other) and
            all(x.exactly_equal(y) for (x, y) in zip(self, other)))

    @classmethod
    def _merge_metadata_dictionaries(cls, dictionaries):
        """
        Helper function for combining variant collections: given multiple
        dictionaries mapping:
             source name -> (variant -> (attribute -> value))

        Returns dictionary with union of all variants and sources.
        """
        # three levels of nested dictionaries!
        #   {source name: {variant: {attribute: value}}}
        combined_dictionary = {}
        for source_to_metadata_dict in dictionaries:
            for source_name, variant_to_metadata_dict in source_to_metadata_dict.items():
                combined_dictionary.setdefault(source_name, {})
                combined_source_dict = combined_dictionary[source_name]
                for variant, metadata_dict in variant_to_metadata_dict.items():
                    combined_source_dict.setdefault(variant, {})
                    combined_source_dict[variant].update(metadata_dict)
        return combined_dictionary

    @classmethod
    def _combine_variant_collections(cls, combine_fn, variant_collections, kwargs):
        """
        Create a single VariantCollection from multiple different collections.

        Parameters
        ----------

        cls : class
            Should be VariantCollection

        combine_fn : function
            Function which takes any number of sets of variants and returns
            some combination of them (typically union or intersection).

        variant_collections : tuple of VariantCollection

        kwargs : dict
            Optional dictionary of keyword arguments to pass to the initializer
            for VariantCollection.
        """
        kwargs["variants"] = combine_fn(*[set(vc) for vc in variant_collections])
        kwargs["source_to_metadata_dict"] = cls._merge_metadata_dictionaries(
            [vc.source_to_metadata_dict for vc in variant_collections])
        kwargs["sources"] = set.union(*([vc.sources for vc in variant_collections]))
        for key, value in variant_collections[0].to_dict().items():
            # If some optional parameter isn't explicitly specified as an
            # argument to union() or intersection() then use the same value
            # as the first VariantCollection.
            #
            # I'm doing this so that the meaning of VariantCollection.union
            # and VariantCollection.intersection with a single argument is
            # the identity function (rather than setting optional parameters
            # to their default values.
            if key not in kwargs:
                kwargs[key] = value
        return cls(**kwargs)

    def union(self, *others, **kwargs):
        """
        Returns the union of variants in a several VariantCollection objects.
        """
        return self._combine_variant_collections(
            combine_fn=set.union,
            variant_collections=(self,) + others,
            kwargs=kwargs)

    def intersection(self, *others, **kwargs):
        """
        Returns the intersection of variants in several VariantCollection objects.
        """
        return self._combine_variant_collections(
            combine_fn=set.intersection,
            variant_collections=(self,) + others,
            kwargs=kwargs)

    _DATAFRAME_COLUMNS = (
        "chr", "start", "ref", "alt", "gene_name", "gene_id",
    )

    def to_dataframe(self):
        """Build a DataFrame from this variant collection."""
        def row_from_variant(variant):
            return OrderedDict([
                ("chr", variant.contig),
                ("start", variant.original_start),
                ("ref", variant.original_ref),
                ("alt", variant.original_alt),
                ("gene_name", ";".join(variant.gene_names)),
                ("gene_id", ";".join(variant.gene_ids))
            ])
        rows = [row_from_variant(v) for v in self]
        # Always return a DataFrame with the expected columns, even
        # when empty, so downstream code (CSV round-trip, joins) doesn't
        # have to special-case len == 0.
        return pd.DataFrame.from_records(rows, columns=self._DATAFRAME_COLUMNS)

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
            explicitly. Pass ``include_header=False`` for legacy
            consumers that don't tolerate comment lines.
        """
        df = self.to_dataframe()
        if not include_header:
            df.to_csv(path, index=False)
            return

        metadata = OrderedDict()
        metadata["varcode_version"] = _varcode_version
        metadata["reference_name"] = self._serialized_reference_name()
        with open(path, "w") as f:
            write_metadata_header(f, metadata)
            df.to_csv(f, index=False)

    def _serialized_reference_name(self):
        """Pick a reference-name string to record in the CSV header.

        Uses the first variant's reference_name. Returns None if the
        collection is empty or has no consistent reference.
        """
        names = {v.reference_name for v in self if v.reference_name}
        if len(names) == 1:
            return next(iter(names))
        return None

    @classmethod
    def from_csv(
            cls,
            path,
            genome=None,
            distinct=True,
            sort_key=variant_ascending_position_sort_key):
        """Rebuild a VariantCollection from a CSV previously written by
        ``VariantCollection.to_csv()``.

        Parameters
        ----------
        path : str
            Path to the CSV file. Lines starting with '#' are treated as
            comments and parsed as ``key=value`` metadata.

        genome : pyensembl.Genome, str, int, or None
            Reference genome to associate with the loaded variants. If
            ``None``, the reference is read from the CSV's metadata
            header (``# reference_name=...``). If neither is available,
            raises ``ValueError``.

        distinct : bool
            Drop duplicate variants (same as the constructor).

        sort_key : callable
            Sort key for the resulting collection.

        Returns
        -------
        VariantCollection
        """
        header = read_metadata_header(path)
        if genome is None:
            genome = header.get("reference_name")
        if genome is None:
            raise ValueError(
                "from_csv needs a reference genome: pass the `genome` "
                "argument explicitly, or write the CSV with "
                "`to_csv(include_header=True)` so `# reference_name=...` "
                "is recorded in the header. Neither was found at %s." % path)

        # Accept either "chr" or "contig" as the contig column so
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
        required = {"start", "ref", "alt"}
        missing = required - set(df.columns)
        if missing:
            raise ValueError(
                "CSV at %s is missing required columns: %s" % (
                    path, sorted(missing)))

        # Extract the columns we need by name and coerce types up front,
        # then iterate with zip. This is robust against column reordering
        # and extra columns (unlike itertuples, which depends on
        # attribute access to valid-identifier column names in fixed
        # positions) and avoids the per-row overhead of iterrows.
        contigs = df[contig_col].astype(str)
        starts = df["start"].astype(int)
        refs = df["ref"].fillna("")
        alts = df["alt"].fillna("")

        variants = [
            Variant(
                contig=contig,
                start=start,
                ref=ref,
                alt=alt,
                genome=genome,
            )
            for contig, start, ref, alt in zip(contigs, starts, refs, alts)
        ]
        return cls(
            variants=variants,
            distinct=distinct,
            sort_key=sort_key,
        )

    def clone_without_ucsc_data(self):
        variants = [v.clone_without_ucsc_data() for v in self]
        return self.clone_with_new_elements(variants)