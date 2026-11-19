#!/bin/bash
# GOW17 uniform time series (tst/inputs/GOW17_uniform_time_series.athinput, 4^3 cells,
# 30 Myr) with Kokkos BDF and the semi-implicit solver. Outputs go to out_<tag>/.
# ATHENA: binary built from ~/Projects/athenak-semi-implicit-base (branch
# semi-implicit-base, uncommitted working tree of 2026-11-11, build-cpu, g++-16).
# Figure: python plot_semi_implicit_vs_bdf.py \
#   --data-dir ~/Documents/Athenak-chem-tests/gow17-semi-implicit-vs-bdf --gpu "Mac CPU" --out gow17_semi_implicit_vs_bdf.png
set -e
ATHENA=${ATHENA:-$HOME/Projects/athenak-semi-implicit-base/build-cpu/src/athena}
IN=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-semi-implicit-vs-bdf}  # run output, outside Dropbox
mkdir -p "$OUT" && cd "$OUT"
run() {  # tag input args...
  local tag=$1 in=$2; shift 2
  rm -rf out_$tag && mkdir -p out_$tag && cd out_$tag
  $ATHENA -i $IN/$in job/basename=$tag "$@" > log 2>&1
  cd ..
}
run bdf_tight  in_bdf.athinput chemistry/kokkos_BDF_rtol=1e-7 chemistry/kokkos_BDF_atol=1e-16
run kokkos_BDF in_bdf.athinput
run q0         in_semi_implicit.athinput
run q1e-1      in_semi_implicit.athinput chemistry/semi_implicit_rtol=0.1 chemistry/semi_implicit_nsub_max=100000
run b9_1e-1    in_semi_implicit.athinput chemistry/semi_implicit_rtol=0.1 chemistry/semi_implicit_yfloor=1e-9 chemistry/semi_implicit_nsub_max=100000
