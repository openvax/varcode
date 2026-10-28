"""Compatibility API for upstream ORFs; see :mod:`varcode.orf` for all ORFs."""

from .orf import (  # noqa: F401
    ORFAnnotation as UORFAnnotation,
    _ORFTranscript,
    _bind,
    _annotate_orf_variants,
    annotate_orf_product,
    predict_orf_effects,
)


def load_uorf_catalog(*args, **kwargs):
    """Load the legacy upstream-ORF catalogue through ORFdata.

    Arguments are forwarded to ``orfdata.load_uorf_catalog``. Requires the
    ``varcode[orf]`` extra. Use :func:`varcode.load_orf_catalog` for the
    general ORF catalogue; this compatibility function keeps its upstream-only
    selection.
    """
    from orfdata import load_uorf_catalog as load
    return load(*args, **kwargs)


class ORFTranscript(_ORFTranscript):
    """Legacy engine view; use ``record.bind(transcript)`` for the base object."""

    def __init__(self, record, transcript, *, genome_build, reference_context=None):
        super().__init__(_bind(record, transcript, genome_build, reference_context=reference_context))


def annotate_uorf_variants(variants, record, transcript, *, genome_build, **kwargs):
    """Annotate variants on a legacy record bound to a supplied transcript.

    Returns an :class:`varcode.ORFAnnotation`. ``genome_build`` identifies the
    record's placement; additional keywords go to the shared ORF annotation
    engine. New callers can bind the reference and use
    :func:`varcode.annotate_orf_variants` directly.
    """
    return _annotate_orf_variants(variants, record, transcript,
                                 genome_build=genome_build, **kwargs)


def annotate_uorf_product(record, transcript, product, *, genome_build, reference_context=None):
    """Compare an observed product with a legacy record/transcript placement.

    Returns a JSON-compatible dictionary. The record is bound using
    ``genome_build`` and optional ``reference_context`` before delegating to
    :func:`varcode.annotate_orf_product`.
    """
    return annotate_orf_product(_bind(record, transcript, genome_build,
                                     reference_context=reference_context), product)


def predict_uorf_effects(variants, *, catalog=None, **kwargs):
    """Return a tuple of ORFAnnotation objects from the upstream-ORF catalogue.

    Omitted ``catalog`` loads :func:`load_uorf_catalog`. Remaining arguments
    are forwarded to :func:`varcode.predict_orf_effects`; supply a catalogue
    explicitly to control its source and reference context.
    """
    return predict_orf_effects(variants,
                              catalog=load_uorf_catalog() if catalog is None else catalog,
                              **kwargs)
