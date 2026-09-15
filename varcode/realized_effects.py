"""Compatibility aliases for the experimental :mod:`varcode.transcript_model`.

New callers should use ``TranscriptModelEffectAnnotator`` and
``predict_transcript_model_effect``. The old imports remain valid.
"""

from .transcript_model import (
    TranscriptModelEffectAnnotator as RealizedEffectAnnotator,
    predict_transcript_model_effect as predict_realized_effect,
)

__all__ = ["RealizedEffectAnnotator", "predict_realized_effect"]
