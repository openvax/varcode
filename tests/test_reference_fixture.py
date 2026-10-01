"""Reference selection must not provision data or depend on the user cache."""

import os
from pathlib import Path
import subprocess
import sys

import pytest


@pytest.mark.parametrize("cache_state", ["empty", "partial", "stale_index"])
def test_grch38_fixture_resolution_never_downloads_or_changes_cache(tmp_path, cache_state):
    # A fresh process avoids PyEnsembl/Varcode memoization and reads the cache
    # environment before either library is imported. The only intercepted
    # operation is the source-acquisition hook implicated in #493.
    script = r'''
import os
from pathlib import Path
from unittest.mock import patch

from pyensembl import EnsemblRelease, Genome, cached_release

cache = Path(os.environ["PYENSEMBL_CACHE_DIR"])
state = os.environ["VARCODE_TEST_CACHE_STATE"]
if state != "empty":
    newer = EnsemblRelease(95)
    gtf = Path(newer.download_cache.cached_path(newer.gtf_url))
    gtf.parent.mkdir(parents=True, exist_ok=True)
    gtf.write_bytes(b"incomplete GTF download")
    if state == "stale_index":
        # A database path alone does not establish a usable installation.
        gtf.with_suffix(".db").write_bytes(b"incomplete database")

def snapshot():
    return {str(p.relative_to(cache)): p.read_bytes() if p.is_file() else None
            for p in cache.rglob("*")}

before = snapshot()
attempts = []

def forbid_source_acquisition(*args, **kwargs):
    attempts.append((args, kwargs))
    raise AssertionError("reference resolution attempted source acquisition")

with patch.object(Genome, "_set_local_paths", forbid_source_acquisition):
    import tests.conftest
    from varcode.reference import infer_genome

    references = ["GRCh38", "hg38", "B38", "GRCh38.p13",
                  "/references/GRCh38.d1.vd1.fa"]
    selected = [infer_genome(name)[0] for name in references]
    assert not attempts, "reference selection tried to acquire missing sources"
    assert all(genome.release == 81 for genome in selected)
    # Explicit choices bypass the test default, including object identity.
    assert infer_genome("GRCh38:95")[0].release == 95
    assert infer_genome(95)[0].release == 95
    supplied = cached_release(93)
    assert infer_genome(supplied)[0] is supplied

assert snapshot() == before
'''
    result = subprocess.run(
        [sys.executable, "-c", script],
        cwd=Path(__file__).resolve().parents[1],
        env=dict(os.environ, PYENSEMBL_CACHE_DIR=str(tmp_path),
                 VARCODE_TEST_CACHE_STATE=cache_state),
        capture_output=True, text=True, timeout=30)
    assert result.returncode == 0, result.stdout + result.stderr
