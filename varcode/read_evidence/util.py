# Copyright (c) 2015. Mount Sinai School of Medicine
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

from __future__ import absolute_import

def alignment_key(pysam_alignment_record):
    '''
    Return the identifying attributes of a `pysam.AlignedSegment` instance.
    This is necessary since these objects do not support a useful notion of
    equality (they compare on identify by default).
    '''
    return (
        read_key(pysam_alignment_record),
        pysam_alignment_record.query_alignment_start,
        pysam_alignment_record.query_alignment_end,
    )

def read_key(pysam_alignment_record):
    '''
    Given a `pysam.AlignedSegment` instance, return the attributes identifying
    the *read* it comes from (not the alignment). There may be more than one
    alignment for a read, e.g. chimeric and secondary alignments.
    '''
    return (
        pysam_alignment_record.query_name,
        pysam_alignment_record.is_duplicate,
        pysam_alignment_record.is_read1,
        pysam_alignment_record.is_read2,
    )
