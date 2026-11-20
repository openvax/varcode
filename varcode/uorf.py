"""Deprecated compatibility API for upstream ORFs; will be removed in varcode 12.0.

Use :mod:`varcode.orf` and ``record.bind(transcript)`` for new callers.
"""

from .effects.common import warn_deprecated
from .orf import (
    ORFAnnotation,
    ORFTranscriptView,
    annotate_bound_orf_variants,
    annotate_orf_product,
    bind_orf_reference,
    predict_orf_effects,
)


class _UORFAnnotationMeta(type(ORFAnnotation)):
    def __instancecheck__(cls, instance):
        return isinstance(instance, ORFAnnotation)

    def __subclasscheck__(cls, subclass):
        return issubclass(subclass, ORFAnnotation)


class UORFAnnotation(ORFAnnotation, metaclass=_UORFAnnotationMeta):
    """Deprecated alias for :class:`~varcode.ORFAnnotation`; will be removed in varcode 12.0."""

    def __new__(cls, *args, **kwargs):
        warn_deprecated(
            "UORFAnnotation is deprecated and will be removed in varcode 12.0; "
            "use ORFAnnotation")
        return ORFAnnotation(*args, **kwargs)


def load_uorf_catalog(*args, **kwargs):
    """Deprecated wrapper for ``orfdata.load_uorf_catalog``; will be removed in varcode 12.0.

    Arguments are forwarded to ``orfdata.load_uorf_catalog``. Use
    ``orfdata.load_uorf_catalog`` directly for the upstream-only selection or
    :func:`varcode.load_orf_catalog` for the general ORF catalogue.
    """
    warn_deprecated(
        "load_uorf_catalog is deprecated and will be removed in varcode 12.0; "
        "use orfdata.load_uorf_catalog or varcode.load_orf_catalog")
    from orfdata import load_uorf_catalog as load
    return load(*args, **kwargs)


class ORFTranscript(ORFTranscriptView):
    """Deprecated engine view; will be removed in varcode 12.0. Use ``record.bind(transcript)``."""

    def __init__(self, record, transcript, *, genome_build, reference_context=None):
        warn_deprecated(
            "ORFTranscript is deprecated and will be removed in varcode 12.0; "
            "use record.bind(transcript)")
        super().__init__(bind_orf_reference(record, transcript, genome_build, reference_context=reference_context))


def annotate_uorf_variants(variants, record, transcript, *, genome_build, **kwargs):
    """Deprecated record-plus-transcript wrapper; will be removed in varcode 12.0.

    Bind the reference with ``record.bind(transcript)`` and use
    :func:`varcode.annotate_orf_variants`.
    """
    warn_deprecated(
        "annotate_uorf_variants is deprecated and will be removed in varcode 12.0; "
        "bind the record and use annotate_orf_variants")
    return annotate_bound_orf_variants(variants, record, transcript,
                                       genome_build=genome_build, **kwargs)


def annotate_uorf_product(record, transcript, product, *, genome_build, reference_context=None):
    """Deprecated record-plus-transcript product wrapper; will be removed in varcode 12.0.

    Bind the reference with ``record.bind(transcript)`` and use
    :func:`varcode.annotate_orf_rna_product` or :func:`varcode.annotate_orf_product`.
    """
    warn_deprecated(
        "annotate_uorf_product is deprecated and will be removed in varcode 12.0; "
        "bind the record and use annotate_orf_rna_product or annotate_orf_product")
    return annotate_orf_product(bind_orf_reference(record, transcript, genome_build,
                                                   reference_context=reference_context), product)


def predict_uorf_effects(variants, *, catalog=None, **kwargs):
    """Deprecated upstream-ORF catalogue wrapper; will be removed in varcode 12.0.

    Use :func:`varcode.predict_orf_effects` with ``catalog=orfdata.load_uorf_catalog()``.
    """
    warn_deprecated(
        "predict_uorf_effects is deprecated and will be removed in varcode 12.0; "
        "use predict_orf_effects with catalog=orfdata.load_uorf_catalog()")
    if catalog is None:
        from orfdata import load_uorf_catalog as load
        catalog = load()
    return predict_orf_effects(variants, catalog=catalog, **kwargs)


__all__ = [
    "ORFTranscript",
    "UORFAnnotation",
    "annotate_uorf_product",
    "annotate_uorf_variants",
    "load_uorf_catalog",
    "predict_uorf_effects",
]
