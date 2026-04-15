#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"

echo "==> Building single-precision version"
gfortran -std=legacy -Wall -Wextra heisen_multi_full.f -o heisen_multi_full

echo "==> Running single-precision version"
./heisen_multi_full | tee heisen_multi_full_output.txt

echo "==> Building double-precision version"
gfortran -std=legacy -Wall -Wextra heisen_multi_double.f -o heisen_multi_double

echo "==> Running double-precision version"
./heisen_multi_double | tee heisen_multi_double_output.txt

echo "==> Done"
echo "Created:"
echo "  heisen_multi_full"
echo "  heisen_multi_double"
echo "  heisen_multi_full_output.txt"
echo "  heisen_multi_double_output.txt"
