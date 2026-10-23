#!/bin/bash
set -o errexit

"${PYTHON:-python}" -m ruff check varcode/ tests/
"${PYTHON:-python}" -m mypy

echo 'Passes ruff and mypy checks'
