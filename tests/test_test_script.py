"""Verify how test.sh invokes pytest, using a stand-in python on PATH."""

import os
from pathlib import Path
import subprocess

import pytest

TEST_SH = Path(__file__).resolve().parents[1] / "test.sh"

# Records pytest's arguments, and whether the --basetemp directory exists
# while "pytest" runs, then exits with $PYTEST_STATUS. "python -c" is
# test.sh's xdist probe, answered with $PROBE_STATUS.
STAND_IN_PYTHON = """#!/bin/sh
if [ "$1" = "-c" ]; then exit "$PROBE_STATUS"; fi
printf "%s\\n" "$@" > "$PYTHON_ARGS"
previous=""
for argument in "$@"; do
    if [ "$previous" = "--basetemp" ] && [ -d "$argument" ]; then
        : > "$PYTHON_ARGS.basetemp-existed"
    fi
    previous="$argument"
done
exit "$PYTEST_STATUS"
"""


def run_test_script(tmp_path, *arguments, probe_status=1, pytest_status=0, **environment):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    interpreter = bin_dir / "python"
    interpreter.write_text(STAND_IN_PYTHON)
    interpreter.chmod(0o755)
    unrelated_pytest = bin_dir / "pytest"
    unrelated_pytest.write_text('#!/bin/sh\nexit 99\n')
    unrelated_pytest.chmod(0o755)
    temp_root = tmp_path / "tmp"
    temp_root.mkdir()
    args_file = tmp_path / "python-args"
    env = dict(os.environ, PATH=str(bin_dir) + os.pathsep + os.environ["PATH"],
               PROBE_STATUS=str(probe_status), PYTEST_STATUS=str(pytest_status),
               PYTHON_ARGS=str(args_file), TMPDIR=str(temp_root),
               TEST_SH_MAX="1", TEST_SH_MIN="1", PER_WORKER_GB="1.5")
    env.pop("PYTEST_ADDOPTS", None)
    env.update(environment)
    result = subprocess.run(["bash", str(TEST_SH), *arguments],
                            env=env, capture_output=True, text=True)
    return result, args_file.read_text().splitlines(), temp_root


@pytest.mark.parametrize("has_xdist", [False, True])
def test_test_script_uses_one_interpreter(tmp_path, has_xdist):
    result, args, _ = run_test_script(
        tmp_path, "-k", "a or b", probe_status=0 if has_xdist else 1, pytest_status=17)
    assert result.returncode == 17, result.stderr
    assert args[:2] == ["-m", "pytest"]
    assert ("-n" in args) is has_xdist
    assert args[-2:] == ["-k", "a or b"]


@pytest.mark.parametrize("pytest_status", [0, 1])
def test_each_run_gets_its_own_temp_directory(tmp_path, pytest_status):
    # pytest's shared $TMPDIR/pytest-of-$USER keeps only three runs, so other
    # runs could delete this one's tmp_path directories (#718).
    result, args, temp_root = run_test_script(tmp_path, pytest_status=pytest_status)
    assert result.returncode == pytest_status, result.stderr
    assert args.count("--basetemp") == 1
    basetemp = Path(args[args.index("--basetemp") + 1])
    assert basetemp.parent == temp_root
    assert basetemp.name.startswith("varcode-pytest.")
    assert (tmp_path / "python-args.basetemp-existed").exists()
    # Removed after passing tests; kept, with its path logged, after failures.
    assert basetemp.exists() is (pytest_status != 0)
    kept = "kept this run's temporary files in %s\n" % basetemp
    assert (kept in result.stderr) is (pytest_status != 0)


@pytest.mark.parametrize("arguments,environment", [
    (["--basetemp", "CHOSEN"], {}),
    (["--basetemp=CHOSEN"], {}),
    ([], {"PYTEST_ADDOPTS": "-q --basetemp=CHOSEN"}),
])
def test_a_chosen_basetemp_is_used_instead(tmp_path, arguments, environment):
    chosen = str(tmp_path / "chosen")
    arguments = [a.replace("CHOSEN", chosen) for a in arguments]
    environment = {k: v.replace("CHOSEN", chosen) for k, v in environment.items()}
    result, args, temp_root = run_test_script(tmp_path, *arguments, **environment)
    assert result.returncode == 0, result.stderr
    assert [a for a in args if a.startswith("--basetemp")] == [
        a for a in arguments if a.startswith("--basetemp")]
    assert list(temp_root.iterdir()) == []
