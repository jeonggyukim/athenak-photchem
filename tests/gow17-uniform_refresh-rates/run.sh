#!/bin/bash
# semi_implicit_refresh_rates true vs false, GOW17 uniform time series (4^3 cells,
# n_H = 109.21 cm^-3, G0 = 1e-6, xi_CR = 2e-16 s^-1, 30 Myr), PR 1 defaults
# otherwise (nsub_max 100), against converged Kokkos BDF (rtol 1e-7, atol 1e-16)
# from the same binary.
# athena_a083ca32: ~/Projects/athenak-semi-implicit-base at a083ca32 (PR 1),
# g++-16, Release, Kokkos Serial, Mac CPU.
set -e
ATHENA=${ATHENA:-$HOME/Documents/Athenak-chem-tests/gow17-refresh-rates/athena_a083ca32}
IN=$(cd "$(dirname "$0")" && pwd)
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-refresh-rates}  # run output, outside Dropbox
mkdir -p "$OUT" && cd "$OUT"
run() {  # tag input args...
  local tag=$1 in=$2; shift 2
  rm -rf out_$tag && mkdir -p out_$tag && cd out_$tag
  /usr/bin/time -p $ATHENA -i $IN/$in job/basename=$tag "$@" > log 2> time.txt
  cd ..
}
run ref in_bdf.athinput chemistry/kokkos_BDF_rtol=1e-7 chemistry/kokkos_BDF_atol=1e-16
for f in 0.3 0.1 0.03; do
  for r in false true; do
    run f${f}_refresh_${r} in_si.athinput chemistry/semi_implicit_dt_sub_frac=$f \
        chemistry/semi_implicit_refresh_rates=$r
  done
done
