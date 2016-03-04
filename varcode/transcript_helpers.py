# Copyright (c) 2016. Mount Sinai School of Medicine
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

from __future__ import print_function, division, absolute_import


def interval_offset_on_transcript(start, end, transcript):
    """
    Given an interval [start:end] and a particular transcript,
    return the start offset of the interval relative to the
    chromosomal positions of the transcript.
    """
    # ensure that start_pos:end_pos overlap with transcript positions
    if start > end:
        raise ValueError(
            "start_pos %d shouldn't be greater than end_pos %d" % (
                start, end))
    if start > transcript.end:
        raise ValueError(
            "Range %d:%d starts after transcript %s (%d:%d)" % (
                start,
                end,
                transcript,
                transcript.start,
                transcript.end))
    if end < transcript.start:
        raise ValueError(
            "Range %d:%d ends before transcript %s (%d:%d)" % (
                start,
                end,
                transcript,
                transcript.start,
                transcript.end))
    # trim the start position to the beginning of the transcript
    if start < transcript.start:
        start = transcript.start
    # trim the end position to the end of the transcript
    if end > transcript.end:
        end = transcript.end
    # return earliest offset into the spliced transcript
    return min(
        transcript.spliced_offset(start),
        transcript.spliced_offset(end))
