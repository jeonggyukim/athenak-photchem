#!/bin/bash
# CO/HCO+ block on and off in the semi-implicit solver, GOW17 uniform time series
# (4^3 cells, n_H = 109.21 cm^-3, G0 = 1e-6, xi_CR = 2e-16 s^-1, 30 Myr), against a
# converged Kokkos BDF reference from the same binary.
# ATHENA: ~/Projects/athenak-chem branch semi-implicit 446ad343, build-cpu-g16
# (g++-16, Release, Serial); that branch has semi_implicit_co_block.
# Figure: python plot_co_hcop_block.py \
#   --data-dir ~/Documents/Athenak-chem-tests/gow17-co-hcop-block --out gow17_co_hcop_block.png
set -e
ATHENA=${ATHENA:-$HOME/Projects/athenak-chem/build-cpu-g16/src/athena}
IN=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-co-hcop-block}  # run output, outside Dropbox
mkdir -p "$OUT" && cd "$OUT"
run() {  # tag input args...
  local tag=$1 in=$2; shift 2
  rm -rf out_$tag && mkdir -p out_$tag && cd out_$tag
  $ATHENA -i $IN/$in job/basename=$tag "$@" > log 2>&1
  cd ..
}
run ref in_bdf.athinput
for f in 0.3 0.1 0.03 0.01; do
  for b in true false; do
    run f${f}_block_${b} in_si.athinput chemistry/semi_implicit_cfl=$f chemistry/semi_implicit_co_block=$b
  done
done
