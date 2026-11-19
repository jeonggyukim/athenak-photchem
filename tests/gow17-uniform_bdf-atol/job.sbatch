#!/usr/bin/env bash
# GOW17 uniform time series: Kokkos BDF at rtol 1e-3 with atol 1e-6 ... 1e-14,
# a converged BDF reference and the semi-implicit solver.
#  - accuracy: 4^3, 0-30 Myr (every cell identical, so 4^3 = one cell)
#  - cost: 64^3 (box scaled with nx, so dx and dt match 4^3), 40 cycles from
#    t = 0 (window A) and 40 cycles from a t = 10 Myr restart (window B)
# Source: ~/athenak-semi-implicit-base (branch semi-implicit-base, 9731520a).
#SBATCH --job-name=bdf-atol
#SBATCH --partition=h200
#SBATCH --gres=gpu:H200:1
#SBATCH --mem=32G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=/gpfs/jeonggyukim/athenak-chem/logs/%x-%j.out

set -eo pipefail
ROOT=$HOME/athenak-semi-implicit-base
D=/gpfs/jeonggyukim/athenak-chem/runs/bdf-atol
IN=$HOME/athenak-bdf-atol-inputs
# Batch jobs start without the module function; load Lmod directly
export MODULEPATH=${MODULEPATH:-/opt/ohpc/pub/modulefiles}
. /opt/ohpc/admin/lmod/lmod/init/bash >/dev/null
module purge
module load cuda/12.6.3
export PATH=$HOME/libs/cmake/bin:$PATH
export OMP_NUM_THREADS=1
nvidia-smi --query-gpu=name --format=csv,noheader

$HOME/ai-notes/scripts/prepare_athenak.py --root $ROOT --device gpu --arch HOPPER90 -j 8
BIN=$ROOT/build-gpu-HOPPER90/src/athena
git -C $ROOT log --oneline -1

# tag solver rtol atol
CONFIGS="bdf_a6:kokkos_BDF:1e-3:1e-6 bdf_a8:kokkos_BDF:1e-3:1e-8 bdf_a10:kokkos_BDF:1e-3:1e-10
bdf_a12:kokkos_BDF:1e-3:1e-12 bdf_a14:kokkos_BDF:1e-3:1e-14 si:semi_implicit:1e-3:1e-6"

run() {  # dir input args...
  local dir=$1 in=$2; shift 2
  rm -rf $dir && mkdir -p $dir && cd $dir
  local t0=$(date +%s.%N)
  $BIN -i $in "$@" > log 2>&1
  echo "$dir wall_s=$(echo "$(date +%s.%N) - $t0" | bc)" | tee -a $D/timing.txt
}

mkdir -p $D && : > $D/timing.txt

# accuracy, 4^3, 0-30 Myr; layout read by plot_semi_implicit_vs_bdf.py:
# acc/out_<tag>/tab/<tag>.hydro_w.*.tab
run $D/acc/out_bdf_ref $IN/in_4.athinput job/basename=bdf_ref chemistry/kokkos_BDF_rtol=1e-7 chemistry/kokkos_BDF_atol=1e-16
for c in $CONFIGS; do IFS=: read tag s r a <<< "$c"
  run $D/acc/out_$tag $IN/in_4.athinput job/basename=$tag chemistry/ode_solver=$s chemistry/kokkos_BDF_rtol=$r chemistry/kokkos_BDF_atol=$a
done

# cost, 64^3, window A: cycles 0-40; logs also as costA/log_<tag>.txt
mkdir -p $D/costA $D/costB
for c in $CONFIGS; do IFS=: read tag s r a <<< "$c"
  run $D/n64A_$tag $IN/in_64.athinput chemistry/ode_solver=$s chemistry/kokkos_BDF_rtol=$r chemistry/kokkos_BDF_atol=$a time/nlim=40
  cp $D/n64A_$tag/log $D/costA/log_$tag.txt
done

# cost, 64^3, window B: 40 cycles from t = 10 Myr
run $D/n64_warm $IN/in_64_warm.athinput chemistry/ode_solver=semi_implicit time/tlim=10
RST=$(ls $D/n64_warm/rst/*.rst | tail -1)
NCYC=$(grep -o "cycle=[0-9]*" $D/n64_warm/log | tail -1 | cut -d= -f2)
echo "restart $RST at cycle $NCYC"
for c in $CONFIGS; do IFS=: read tag s r a <<< "$c"
  dir=$D/n64B_$tag; rm -rf $dir && mkdir -p $dir && cd $dir
  t0=$(date +%s.%N)
  $BIN -r $RST chemistry/ode_solver=$s chemistry/kokkos_BDF_rtol=$r chemistry/kokkos_BDF_atol=$a time/tlim=30 time/nlim=$((NCYC + 40)) > log 2>&1
  echo "$dir wall_s=$(echo "$(date +%s.%N) - $t0" | bc)" | tee -a $D/timing.txt
  cp $dir/log $D/costB/log_$tag.txt
done
echo done
