#!/bin/bash
# Cost and accuracy of refresh_rates with the substep cap lifted.
# (1) accuracy, 4^3, nsub_max = 100000: out_nocap_f*_refresh_*
# (2) cost, 16^3 with the same cell size (box 0.25 pc), tab output off, hst on:
#     out_cost_f*_refresh_*; nsub_max 100 (PR 1 default) and 100000.
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
for f in 0.3 0.1 0.03; do
  for r in false true; do
    run nocap_f${f}_refresh_${r} in_si.athinput chemistry/semi_implicit_dt_sub_frac=$f \
        chemistry/semi_implicit_refresh_rates=$r chemistry/semi_implicit_nsub_max=100000
    for n in 100 100000; do
      run cost_n${n}_f${f}_refresh_${r} in_si_cost.athinput \
          chemistry/semi_implicit_dt_sub_frac=$f chemistry/semi_implicit_refresh_rates=$r \
          chemistry/semi_implicit_nsub_max=$n
    done
  done
done
