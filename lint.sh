#!/bin/bash
set -o errexit

ruff check varcode tests
echo 'Passes ruff check'
