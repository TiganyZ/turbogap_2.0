#!/usr/bin/env bash
# Smoke-tests the scenarios under tests_manual/ against the built turbogap
# binary. For each subdirectory containing an `input` file, runs the binary
# in both serial and 2-rank MPI form, in an isolated scratch copy, under a
# wall-clock timeout.
#
# There is no golden/reference output to diff against here - MC/MD runs are
# stochastic (mostly unseeded) and the scenarios are long-running physics
# simulations, not unit tests. So "pass" means: the run either finishes
# cleanly (exit 0) or is still going with no error when the timeout fires,
# and nothing in its output looks like a crash (Fortran runtime error,
# segfault, MPI job abort, ...). This catches build breaks, crashes, and
# regressions like out-of-bounds array accesses; it does not catch silent
# numerical wrongness.
#
# Usage: tests_manual/run_manual_tests.sh [timeout_seconds]
# Env:   TURBOGAP_TEST_NP (default "1 2") - space-separated rank counts to try
set -u

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TURBOGAP="$ROOT/bin/turbogap"
TIMEOUT="${1:-60}"
RANKS="${TURBOGAP_TEST_NP:-1 2}"
WORKDIR="$(mktemp -d "${TMPDIR:-/tmp}/turbogap_test.XXXXXX")"
trap 'rm -rf "$WORKDIR"' EXIT

if [ ! -x "$TURBOGAP" ]; then
   echo "error: $TURBOGAP not found or not executable - run 'make' first" >&2
   exit 1
fi

# Files a run itself writes out - never symlink/keep these from the source
# tree, so a test run can't clobber committed reference data in tests_manual/.
OUTPUT_GLOBS=(trajectory_out.xyz thermo.log mc.log mc_all.xyz '*.log')

n_pass=0
n_fail=0
declare -a fail_labels=()

run_one() {
   local test_dir="$1" mode="$2" np="$3"
   local rel_name label run_dir log status pattern

   rel_name="${test_dir#"$ROOT"/tests_manual/}"
   label="$rel_name [$mode, np=$np]"
   run_dir="$WORKDIR/$(echo "$rel_name.$mode.$np" | tr '/ ' '__')"

   mkdir -p "$run_dir"
   cp -r "$test_dir"/. "$run_dir"/
   for pattern in "${OUTPUT_GLOBS[@]}"; do
      rm -f "$run_dir"/$pattern
   done

   log="$run_dir/run.log"
   (
      cd "$run_dir" || exit 1
      if [ "$np" -gt 1 ]; then
         timeout "$TIMEOUT" mpirun -np "$np" "$TURBOGAP" "$mode"
      else
         timeout "$TIMEOUT" "$TURBOGAP" "$mode"
      fi
   ) >"$log" 2>&1
   status=$?

   if grep -qE "Fortran runtime error|Error termination|SIGSEGV|core dumped|forrtl|exited with non-zero status" "$log"; then
      echo "FAIL  $label (runtime error - see $log)"
      n_fail=$((n_fail + 1))
      fail_labels+=("$label")
      return
   fi

   # 0 = finished cleanly; 124 = our own timeout fired while it was still
   # running with no error seen - a smoke-test pass, not a hang/crash.
   if [ "$status" -eq 0 ] || [ "$status" -eq 124 ]; then
      echo "PASS  $label"
      n_pass=$((n_pass + 1))
   else
      echo "FAIL  $label (exit code $status - see $log)"
      n_fail=$((n_fail + 1))
      fail_labels+=("$label")
   fi
}

# tests_manual/<category>/<scenario>/input - category names the turbogap mode.
declare -A MODE_FOR_CATEGORY=(
   [test_md]=md
   [test_mc]=mc
   [test_predict]=predict
)

shopt -s nullglob
input_files=("$ROOT"/tests_manual/*/*/input "$ROOT"/tests_manual/*/input)
if [ "${#input_files[@]}" -eq 0 ]; then
   echo "No tests_manual/*/*/input files found - nothing to run."
   exit 0
fi

for input_file in "${input_files[@]}"; do
   test_dir="$(dirname "$input_file")"
   rel_name="${test_dir#"$ROOT"/tests_manual/}"
   # tests_manual/<category>/<scenario>/input -> category is the first path
   # component; tests_manual/<category>/input -> category is rel_name itself.
   case "$rel_name" in
      */*) category="${rel_name%%/*}" ;;
      *) category="$rel_name" ;;
   esac
   mode="${MODE_FOR_CATEGORY[$category]:-}"
   if [ -z "$mode" ]; then
      echo "SKIP  ${test_dir#"$ROOT"/tests_manual/} (unknown turbogap mode for category '$category')"
      continue
   fi
   for np in $RANKS; do
      run_one "$test_dir" "$mode" "$np"
   done
done

echo
echo "$n_pass passed, $n_fail failed"
if [ "$n_fail" -gt 0 ]; then
   echo "Failed: ${fail_labels[*]}"
   exit 1
fi
