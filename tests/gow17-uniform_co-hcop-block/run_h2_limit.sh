#!/bin/bash
# As run.sh, but the substep is limited by H2 and energy only (CO removed from the
# step-limiting species: n_step_limiting_species = 1 in gow17.hpp of a scratch
# worktree of semi-implicit 446ad343, branch scratch-h2-only-step-limit, not
# committed). Reference: out_ref from run.sh.
set -e
ATHENA=${ATHENA:?set ATHENA to the H2-only build}
IN=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-co-hcop-block}  # run output, outside Dropbox
mkdir -p "$OUT" && cd "$OUT"
run() {  # tag input args...
  local tag=$1 in=$2; shift 2
  rm -rf out_$tag && mkdir -p out_$tag && cd out_$tag
  $ATHENA -i $IN/$in job/basename=$tag "$@" > log 2>&1
  cd ..
}
for f in 0.3 0.1 0.03 0.01; do
  for b in true false; do
    run h2lim_f${f}_block_${b} in_si.athinput chemistry/semi_implicit_cfl=$f chemistry/semi_implicit_co_block=$b
  done
done
