#!/bin/bash
# Decaying turbulence with GOW17 chemistry: AthenaK (semi-implicit, Kokkos BDF)
# against Athena++ CVODE, all started from one shared VTK snapshot.
#
#   bash tests/turb-box-athenak-vs-athenapp/run.sh                 every stage, in order
#   bash tests/turb-box-athenak-vs-athenapp/run.sh si cvode join   only the named stages
#   SMOKE=1 bash tests/turb-box-athenak-vs-athenapp/run.sh         every stage, 20 cycles each
#
# Stages:
#   drive  AthenaK driven turbulence, hydro only, t = 0 to 12 (1 rank, ~5 min at 64^3)
#   ic     last driving dump -> shared VTK snapshot (and its n_H = 1 copy for CASE=n1)
#   si     AthenaK, semi-implicit chemistry, t = 0 to 2.0454 (NP ranks, ~30 s)
#   bdf    AthenaK, Kokkos BDF at the CVODE tolerances, same interval (NP ranks)
#   cvode  Athena++, CVODE chemistry, same interval (NP ranks, ~18 min)
#   join   merge the Athena++ per-MeshBlock VTK files of each dump into joined/
#   plot   history from the .hst files, history from the dumps, and slice
#          comparisons (final frame and movie)
# Times are for 8 ranks on an Apple-silicon MacBook Pro.
#
# Settings, from the environment:
#   NX    resolution, NX^3 cells (default 64). The driving run is one MeshBlock on
#         one rank; the chemistry runs use MB^3 MeshBlocks (default min(32, NX/2)).
#   CASE  n10 (default): mean n_H = 10 cm^-3, G0 = 1e-6;
#         n1: the same snapshot with density and pressure x 0.1 (n_H = 1 cm^-3),
#         G0 = 1, floors x 0.1.
# Paths and rank count come from machines/<MACHINE>.sh. MACHINE defaults to the
# author's profile on the author's hosts (mac, grammar, syntax) and to default
# everywhere else; CHEM_MACHINE names any other profile file. Every setting can
# be overridden from the environment (ATHENAK_ROOT, ATHENAPP_ROOT, AKBIN, APBIN,
# CHEMDATA, PYTHON, NP).
#
# Output goes to $CHEMDATA (for NX other than 64, <data>-nx<NX> unless CHEMDATA is
# set). A stage refuses to overwrite a run directory that already holds run.out.
# Under Slurm, submit from the repository root: Slurm runs a spooled copy of this
# file, which cannot locate its own directory.
set -e
NX=${NX:-64}
MB=${MB:-$(( NX / 2 < 32 ? NX / 2 : 32 ))}
CASE=${CASE:-n10}

if [ -n "$SLURM_JOB_ID" ]; then
  HERE=$SLURM_SUBMIT_DIR/tests/turb-box-athenak-vs-athenapp
else
  HERE=$(cd "$(dirname "$0")" && pwd)
fi

if [ -z "$MACHINE" ]; then
  case $(hostname -s) in
    Jeong-Gyus-MacBook-Pro*) MACHINE=mac ;;
    grammar*) MACHINE=grammar ;;
    syntax|syn[0-9]*) MACHINE=syntax ;;
    *) MACHINE=default ;;
  esac
fi
source "${CHEM_MACHINE:-$HERE/machines/$MACHINE.sh}"
[ "$NX" = 64 ] || [ -n "$CHEMDATA" ] || DATA=$DATA-nx$NX
# The python scripts find AthenaK's bin_convert.py through this.
export ATHENAK_ROOT=$AK

# Resolution and case, as command-line overrides of the three inputs
IC=shared_ic$NX.vtk
DRIVE_ARGS=(mesh/nx1=$NX mesh/nx2=$NX mesh/nx3=$NX
            meshblock/nx1=$NX meshblock/nx2=$NX meshblock/nx3=$NX)
case $CASE in
  n10) TAG=turb; CASE_IC=$IC; AK_CASE=(); AP_CASE=() ;;
  n1)  TAG=turb_n1; CASE_IC=shared_ic${NX}_n1.vtk
       AK_CASE=(hydro/dfloor=1e-2 hydro/pfloor=1e-1 chemistry/radiation_G0=1.0)
       AP_CASE=(hydro/dfloor=1e-2 hydro/pfloor=1e-1 chem_radiation/G0=1.0) ;;
  *) echo "ERROR  CASE must be n10 or n1" >&2; exit 1 ;;
esac
CHEM_ARGS=(mesh/nx1=$NX mesh/nx2=$NX mesh/nx3=$NX
           meshblock/nx1=$MB meshblock/nx2=$MB meshblock/nx3=$MB problem/vtkfile=$CASE_IC)

# The smoke test runs every stage for 20 cycles, into its own directory, with
# figures beside the data rather than in the repository.
OVR=()
FIGDIR=${FIGDIR:-$HERE/../../figures}
if [ -n "$SMOKE" ]; then
  [ -n "$CHEMDATA" ] || DATA=$DATA-smoke
  OVR=(time/nlim=20)
  FIGDIR=$DATA/figures
fi

STAGES=${*:-drive ic si bdf cvode join plot}
want() { case " $STAGES " in *" $1 "*) return 0 ;; esac; return 1; }

mkdir -p "$DATA"

# A figure is reproducible only if the code that made its numbers is recorded.
# The working trees may carry uncommitted changes, so the commit alone is not
# enough: save the diff and the untracked sources too.
stamp() {
  local name=$1 repo=$2
  {
    echo "$name: $(git -C "$repo" rev-parse --short HEAD)" \
         "on $(git -C "$repo" rev-parse --abbrev-ref HEAD)"
    git -C "$repo" status --porcelain | sed 's/^/  /'
  } >> "$DATA/versions.txt"
  git -C "$repo" diff HEAD > "$DATA/$name.diff"
  git -C "$repo" ls-files --others --exclude-standard -- src > "$DATA/$name-untracked.txt"
  if [ -s "$DATA/$name-untracked.txt" ]; then
    tar -C "$repo" -cf "$DATA/$name-untracked.tar" -T "$DATA/$name-untracked.txt"
  fi
}

check_fresh() {
  local newer
  newer=$(find "$AK/src" -newer "$AKBIN" \( -name '*.cpp' -o -name '*.hpp' \) | head -1)
  if [ -n "$newer" ]; then
    echo "ERROR  $newer is newer than $AKBIN; rebuild first." >&2
    exit 1
  fi
}

check_ap() {
  local cfg
  cfg=$("$APBIN" -c)
  if ! echo "$cfg" | grep -qE 'Riemann solver: +hllc'; then
    echo "ERROR  $APBIN is not built with --flux=hllc; AthenaK uses hllc." >&2
    exit 1
  fi
  if ! echo "$cfg" | grep -qE 'Problem generator: +read_vtk'; then
    echo "ERROR  $APBIN is not built with --prob=read_vtk." >&2
    exit 1
  fi
  echo "$cfg" > "$DATA/athena-pp-configure.txt"
}

prepare() {  # run directory, input file, initial snapshot to copy in (optional)
  local dir=$DATA/$1
  if [ -e "$dir/run.out" ]; then
    echo "ERROR  $dir already holds a run; move it aside first." >&2
    exit 1
  fi
  mkdir -p "$dir"
  cp -f "$HERE/inputs/$2" "$dir/"
  [ -z "$3" ] || cp -f "$DATA/$3" "$dir/"
}

run() {
  local dir=$DATA/$1 np=$2
  shift 2
  echo "=== $1 in $dir ($np ranks)"
  ( cd "$dir" && { time mpirun -np "$np" "$@" > run.out 2> run.err ; } 2> time.txt )
  tail -3 "$dir/time.txt"
}

echo "date: $(date '+%Y-%m-%d %H:%M'), host: $(hostname -s), machine: $MACHINE," \
  "NX: $NX, CASE: $CASE, stages: $STAGES ${SMOKE:+(smoke)}" >> "$DATA/versions.txt"
stamp athenak "$AK"
stamp athena-pp "$AP"

if want drive; then
  check_fresh
  prepare drive_run drive.athinput
  run drive_run 1 "$AKBIN" -i drive.athinput "${DRIVE_ARGS[@]}" "${OVR[@]}"
fi

if want ic; then
  last=$(ls "$DATA"/drive_run/bin/drive.hydro_w.*.bin | tail -1)
  "$PY" "$HERE/athenak_bin_to_vtk.py" "$last" "$DATA/$IC"
  echo "$IC from $last" >> "$DATA/versions.txt"
  if [ "$CASE" = n1 ]; then
    "$PY" "$HERE/scale_vtk_density.py" "$DATA/$IC" "$DATA/$CASE_IC" --factor 0.1
    echo "$CASE_IC from $IC, density and pressure x 0.1" >> "$DATA/versions.txt"
  fi
fi

if want si; then
  check_fresh
  prepare ${TAG}_si_run turb_athenak.athinput "$CASE_IC"
  run ${TAG}_si_run "$NP" "$AKBIN" -i turb_athenak.athinput job/basename=${TAG}_si \
    chemistry/ode_solver=semi_implicit "${CHEM_ARGS[@]}" "${AK_CASE[@]}" "${OVR[@]}"
fi

if want bdf; then
  check_fresh
  prepare ${TAG}_bdf_run turb_athenak.athinput "$CASE_IC"
  run ${TAG}_bdf_run "$NP" "$AKBIN" -i turb_athenak.athinput job/basename=${TAG}_bdf \
    chemistry/ode_solver=kokkos_BDF "${CHEM_ARGS[@]}" "${AK_CASE[@]}" "${OVR[@]}"
fi

if want cvode; then
  check_ap
  prepare ${TAG}_cvode_run turb_athenapp.athinput "$CASE_IC"
  run ${TAG}_cvode_run "$NP" "$APBIN" -i turb_athenapp.athinput job/problem_id=${TAG}_cvode \
    "${CHEM_ARGS[@]}" "${AP_CASE[@]}" "${OVR[@]}"
fi

if want join; then
  dir=$DATA/${TAG}_cvode_run
  base=${TAG}_cvode
  nblock=$(ls "$dir"/$base.block*.out1.00000.vtk | wc -l)
  mkdir -p "$dir/joined"
  for f in "$dir"/$base.block0.out1.*.vtk; do
    n=${f##*.out1.}
    blocks=()
    for (( b = 0; b < nblock; b++ )); do blocks+=("$dir/$base.block$b.out1.$n"); done
    "$AP/vis/vtk/join_vtk++" -o "$dir/joined/$base.block0.out1.$n" "${blocks[@]}" > /dev/null
  done
  echo "joined $(ls "$dir/joined" | wc -l | tr -d ' ') dumps of $((nblock)) MeshBlocks each"
fi

if want plot; then
  mkdir -p "$FIGDIR"
  "$PY" "$HERE/../../scripts/plot_turb_hst.py" \
    --ak "$DATA/${TAG}_si_run/${TAG}_si.hydro.hst" \
    --ap "$DATA/${TAG}_cvode_run/${TAG}_cvode.hst" \
    --out "$FIGDIR/fig_${TAG}_hst.png"
  "$PY" "$HERE/../../scripts/plot_turb_history_dumps.py" \
    --ak-si-dir "$DATA/${TAG}_si_run" --ak-si-base ${TAG}_si \
    --ak-bdf-dir "$DATA/${TAG}_bdf_run" --ak-bdf-base ${TAG}_bdf \
    --ap-cvode-dir "$DATA/${TAG}_cvode_run" --ap-cvode-base ${TAG}_cvode \
    --out "$FIGDIR/fig_${TAG}_history_dumps.png"
  "$PY" "$HERE/../../scripts/plot_turb_evolution.py" \
    --ak-dir "$DATA/${TAG}_si_run" --ak-base ${TAG}_si \
    --ap-dir "$DATA/${TAG}_cvode_run/joined" --ap-base ${TAG}_cvode \
    --outdir "$FIGDIR"
fi
