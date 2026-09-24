"""Verify test.sh uses the interpreter whose optional plugins it probes."""

import os
from pathlib import Path
import subprocess

import pytest


@pytest.mark.parametrize("has_xdist", [False, True])
def test_test_script_uses_one_interpreter(tmp_path, has_xdist):
    interpreter = tmp_path / "python"
    interpreter.write_text(
        '#!/bin/sh\n'
        'if [ "$1" = "-c" ]; then exit "$PROBE_STATUS"; fi\n'
        'printf "%s\\n" "$@" > "$PYTHON_ARGS"\n'
        'exit 17\n'
    )
    interpreter.chmod(0o755)
    unrelated_pytest = tmp_path / "pytest"
    unrelated_pytest.write_text('#!/bin/sh\nexit 99\n')
    unrelated_pytest.chmod(0o755)
    args_file = tmp_path / "python-args"
    env = dict(os.environ, PATH=str(tmp_path) + os.pathsep + os.environ["PATH"],
               PROBE_STATUS="0" if has_xdist else "1", PYTHON_ARGS=str(args_file),
               TEST_SH_MAX="1", TEST_SH_MIN="1", PER_WORKER_GB="1.5")
    result = subprocess.run(
        ["bash", str(Path(__file__).resolve().parents[1] / "test.sh"), "-k", "a or b"],
        env=env, capture_output=True, text=True)
    assert result.returncode == 17, result.stderr
    args = args_file.read_text().splitlines()
    assert args[:2] == ["-m", "pytest"]
    assert ("-n" in args) is has_xdist
    assert args[-2:] == ["-k", "a or b"]
