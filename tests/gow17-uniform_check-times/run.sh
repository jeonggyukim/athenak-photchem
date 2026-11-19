#!/bin/bash
# GOW17 uniform (4^3): one run per check time, stopped at tlim (hit exactly), for the
# converged Kokkos BDF reference (rtol 1e-7, atol 1e-16), default Kokkos BDF
# (rtol 1e-3, atol 1e-6) and semi_implicit at its defaults (dt_sub_frac 0.1).
# Conditions (name: start, n_H [cm^-3], xi_CR [s^-1], G0; check times [Myr]):
#   cur: all species 1e-6, 109.21, 2e-16, 1e-6; 0.1 0.4 3 30
#   B:   molecular, 1e2, 2e-16, 1;             1e-4 1e-3 0.01 0.1
#   C:   atomic, 1e4, 2e-16, 1e-6;             0.02 0.06 0.15 1
#   D:   atomic, 1e3, 2e-15, 1e-6;             0.5 0.8 1.1 5
# ATHENA: binary of the solver under test; output to OUT (outside Dropbox).
# ONLY=<name> runs one condition. Errors: python check_errors.py (in OUT).
set -e
IN=$(cd "$(dirname "$0")" && pwd)
: "${ATHENA:?set ATHENA}"
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-uniform_check-times}
mkdir -p "$OUT/out" && cd "$OUT/out"
run() {  # tag input tlim args...
  local tag=$1 in=$2 t=$3; shift 3
  rm -rf $tag && mkdir $tag && cd $tag
  $ATHENA -i $IN/$in job/basename=run time/tlim=$t output1/dcycle=1000000 "$@" > log 2>&1
  cd ..
}
case_runs() {  # name start nH cr g0 times...
  local name=$1 st=$2 n=$3 cr=$4 g=$5; shift 5
  local a=(problem/n_H=$n chemistry/radiation_CR=$cr chemistry/radiation_G0=$g)
  for t in "$@"; do
    run ${name}_t${t}_ref bdf_${st}.athinput $t "${a[@]}" \
        chemistry/kokkos_BDF_rtol=1e-7 chemistry/kokkos_BDF_atol=1e-16
    run ${name}_t${t}_bdf bdf_${st}.athinput $t "${a[@]}"
    run ${name}_t${t}_si si_${st}.athinput $t "${a[@]}"
  done
}
[ -n "$ONLY" ] && [ "$ONLY" != cur ] || case_runs cur 1e-6 1.0921e2 2e-16 1e-6 0.1 0.4 3 30
[ -n "$ONLY" ] && [ "$ONLY" != B ] || case_runs B molecular 1e2 2e-16 1 1e-4 1e-3 0.01 0.1
[ -n "$ONLY" ] && [ "$ONLY" != C ] || case_runs C atomic 1e4 2e-16 1e-6 0.02 0.06 0.15 1
[ -n "$ONLY" ] && [ "$ONLY" != D ] || case_runs D atomic 1e3 2e-15 1e-6 0.5 0.8 1.1 5
