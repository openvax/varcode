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

"""
Time how long it takes to open a VCF.

Run as:
    python -m profile -s cumtime %(prog)s

to get profiling output.

"""
import argparse
import time

import varcode

parser = argparse.ArgumentParser(description=__doc__)

parser.add_argument(
    "path", help="Path or URL to VCF")

parser.add_argument(
    "--profile", action="store_true",
    default=False,
    help="Run in a profiler.")

parser.add_argument(
    "--no-info-field",
    dest="info_field",
    action="store_false",
    default=True)

parser.add_argument(
    "--pyvcf",
    help="use pyvcf implementation",
    action="store_true",
    default=False)

def run():
    args = parser.parse_args()

    extra_args = {}
    if not args.info_field:
        extra_args["include_info"] = False

    start = time.time()

    if args.pyvcf:
        result = varcode.load_vcf(
            args.path,
            allow_extended_nucleotides=True)
    else:
        result = varcode.load_vcf_fast(
            args.path,
            allow_extended_nucleotides=True,
            **extra_args)

    print("Loaded %d variants in %0.3f sec. " % (
        len(result), time.time() - start))
    print(result.to_string(limit=5))

if __name__ == '__main__':
    run()
