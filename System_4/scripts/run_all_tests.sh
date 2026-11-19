#!/usr/bin/env bash
# Builds System 4 and runs every test. Exit code != 0 if anything fails.
#   scripts/run_all_tests.sh              build + core-rule check + all unit tests
#   scripts/run_all_tests.sh --analysis   additionally run the 60 s analysis + plots
set -u
cd "$(dirname "$0")/.."
mkdir -p build
fail=0

echo "### build (-std=c99 -Wall -Wextra -Wpedantic)"
if ! make -s all 2>&1 | tee build/make.log | grep -E "warning|error"; then :; fi
if grep -qE "warning|error" build/make.log; then echo "BUILD HAS WARNINGS/ERRORS"; fail=1; fi
[ -x build/bin/run_analysis ] || { echo "build failed"; exit 1; }

echo; echo "### coding-rule check"
python3 scripts/check_core_rules.py || fail=1

for t in test_fixed_math test_divider test_cordic test_adapter test_quant_stress test_models; do
    echo; echo "### $t"
    if ! "build/bin/$t"; then echo "*** $t FAILED"; fail=1; fi
done

echo; echo "### system4_proof (best-configuration dummy-test proof)"
if ! build/bin/system4_proof results; then echo "*** system4_proof FAILED"; fail=1; fi

if [ "${1:-}" = "--analysis" ]; then
    echo; echo "### analysis (60 s scenarios)"
    build/bin/run_analysis results 60 || fail=1
    python3 scripts/plot_results.py results || fail=1
fi

echo
if [ $fail -eq 0 ]; then echo "ALL TESTS PASSED"; else echo "SOME TESTS FAILED"; fi
exit $fail
