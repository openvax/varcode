"""Deprecated compatibility aliases for :mod:`varcode.transcript_model`.

These aliases emit ``DeprecationWarning`` and will be removed in varcode 12.0;
use :class:`~varcode.TranscriptModelEffectAnnotator` and
:func:`~varcode.predict_transcript_model_effect`.
"""

from functools import wraps

from .effects.common import warn_deprecated
from .transcript_model import (
    TranscriptModelEffectAnnotator,
    predict_transcript_model_effect,
)


class RealizedEffectAnnotator(TranscriptModelEffectAnnotator):
    """Deprecated alias for :class:`~varcode.TranscriptModelEffectAnnotator`; will be removed in varcode 12.0."""

    def __init__(self, *args, **kwargs):
        warn_deprecated(
            "RealizedEffectAnnotator is deprecated and will be removed in "
            "varcode 12.0; use TranscriptModelEffectAnnotator")
        super().__init__(*args, **kwargs)

    def __reduce__(self):
        return (TranscriptModelEffectAnnotator, ())


@wraps(predict_transcript_model_effect)
def predict_realized_effect(*args, **kwargs):
    """Deprecated alias for :func:`~varcode.predict_transcript_model_effect`; will be removed in varcode 12.0."""
    warn_deprecated(
        "predict_realized_effect is deprecated and will be removed in "
        "varcode 12.0; use predict_transcript_model_effect")
    return predict_transcript_model_effect(*args, **kwargs)


__all__ = ["RealizedEffectAnnotator", "predict_realized_effect"]
