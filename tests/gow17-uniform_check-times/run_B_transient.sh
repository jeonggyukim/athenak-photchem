#!/bin/bash
# Condition B (molecular start, n_H 1e2, G0 1): the photodissociation transient lies
# inside the first hydro step (0.0097 Myr), where no regular output exists. One
# tlim run per time, 1e-5 to 1e-2 Myr, converged BDF and semi_implicit, for the
# plot only. ATHENA as in run.sh; output to OUT/transient_B.
set -e
IN=$(cd "$(dirname "$0")" && pwd)
: "${ATHENA:?set ATHENA}"
OUT=${OUT:-$HOME/Documents/Athenak-chem-tests/gow17-uniform_check-times}
mkdir -p "$OUT/transient_B" && cd "$OUT/transient_B"
a=(problem/n_H=1e2 chemistry/radiation_CR=2e-16 chemistry/radiation_G0=1 output1/dcycle=1000000)
for t in 1e-5 2e-5 5e-5 1e-4 2e-4 5e-4 1e-3 2e-3 5e-3 1e-2; do
  for s in ref si; do
    d=t${t}_$s; rm -rf $d && mkdir $d && cd $d
    if [ $s = ref ]; then
      $ATHENA -i $IN/bdf_molecular.athinput job/basename=run time/tlim=$t "${a[@]}" \
          chemistry/kokkos_BDF_rtol=1e-7 chemistry/kokkos_BDF_atol=1e-16 > log 2>&1
    else
      $ATHENA -i $IN/si_molecular.athinput job/basename=run time/tlim=$t "${a[@]}" > log 2>&1
    fi
    cd ..
  done
done
