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
    from orfdata import load_uorf_catalog as load
    return load(*args, **kwargs)


class ORFTranscript(_ORFTranscript):
    """Legacy engine view; use ``record.bind(transcript)`` for the base object."""

    def __init__(self, record, transcript, *, genome_build):
        super().__init__(_bind(record, transcript, genome_build))


def annotate_uorf_variants(variants, record, transcript, *, genome_build, **kwargs):
    return _annotate_orf_variants(variants, record, transcript,
                                 genome_build=genome_build, **kwargs)


def annotate_uorf_product(record, transcript, product, *, genome_build):
    return annotate_orf_product(_bind(record, transcript, genome_build), product)


def predict_uorf_effects(variants, *, catalog=None, **kwargs):
    return predict_orf_effects(variants,
                              catalog=load_uorf_catalog() if catalog is None else catalog,
                              **kwargs)
