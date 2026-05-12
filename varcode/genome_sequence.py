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

"""Reference sequence lookup against a pyensembl :class:`Genome`.

varcode's reference is the pyensembl ``Genome`` object the user passes
through ``genome=``. By default that genome ships transcript and
protein FASTAs but not the chromosome FASTA — so any feature that
needs raw bases at arbitrary genomic positions (intronic, intergenic,
flanking) has to go through this module's tiered lookup:

1. **Tier 1** — if the genome has an attached chromosome FASTA via
   :func:`attach_genome_fasta`, read from it.
2. **Tier 2** — fall back to pyensembl transcript cDNA via
   ``transcript.spliced_offset()`` for any transcript covering the
   position. Reverse-complements for ``-`` strand transcripts so the
   result is always on the ``+`` strand.
3. **Tier 3** — return ``""``. The caller decides whether empty is an
   error or just "no shift possible / no candidates / etc."

Internal API: feature code calls :func:`reference_base` or
:func:`reference_range`. User-facing API: :func:`attach_genome_fasta`.

Tracked in openvax/varcode#372.
"""

_VARCODE_FASTA_ATTR = "_varcode_fasta"


def attach_genome_fasta(genome, fasta_or_path):
    """Attach a chromosome FASTA to a pyensembl ``Genome`` so varcode
    can read raw reference bases at arbitrary positions.

    The default ``pyensembl install`` ships transcript and protein
    FASTAs only. Features that need intronic, intergenic, or flanking
    bases (left-alignment, cryptic-exon scoring, sequence-aware
    splice prediction) require this attachment.

    Parameters
    ----------
    genome : pyensembl.Genome
        The genome to extend. Mutated in place — the FASTA is stored
        on a varcode-specific underscore attribute.
    fasta_or_path : str, path-like, or FASTA object
        One of:

        * A path string (``"/path/to/GRCh38.fa"``) — opened with
          pyfaidx internally.
        * A ``pyfaidx.Fasta`` instance — used as-is.
        * Any object supporting ``fa[contig][start:end].seq`` —
          tests, custom storage, mmap'd files, etc.

    Examples
    --------

    >>> import varcode
    >>> import pyfaidx
    >>> from pyensembl import EnsemblRelease
    >>> g = EnsemblRelease(81)
    >>> varcode.attach_genome_fasta(g, "/path/to/GRCh38.fa")  # doctest: +SKIP
    >>> vc = varcode.load_vcf("tumor.vcf", genome=g)          # doctest: +SKIP
    """
    raise NotImplementedError(
        "varcode.attach_genome_fasta — implementation in progress, "
        "tracked in openvax/varcode#372")


def reference_base(genome, contig, position):
    """Return the ``+`` strand reference base at ``(contig, position)``.

    1-based inclusive coordinates matching pyensembl. Returns ``""``
    when no tier of the lookup covers this position.

    See module docstring for the tiered fallback rules.
    """
    raise NotImplementedError(
        "varcode.genome_sequence.reference_base — implementation in "
        "progress, tracked in openvax/varcode#372")


def reference_range(genome, contig, start, end):
    """Return the ``+`` strand reference sequence over ``[start, end]``.

    1-based inclusive coordinates. Returns ``""`` if *any* position in
    the range is not covered by the current tier — the lookup is
    all-or-nothing, so callers see a partial range only when they ask
    for one explicitly.

    See module docstring for the tiered fallback rules.
    """
    raise NotImplementedError(
        "varcode.genome_sequence.reference_range — implementation in "
        "progress, tracked in openvax/varcode#372")
