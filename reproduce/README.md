# Decaying turbulence with GOW17 chemistry: AthenaK vs Athena++

This test runs the same turbulent box with GOW17 chemistry in AthenaK and in
Athena++ and compares the two, cell by cell. AthenaK uses its semi-implicit
solver (or Kokkos BDF); Athena++ uses CVODE.

The two codes drive turbulence differently, so the box is driven once, in
AthenaK, and both codes then start from that same snapshot. From there each
code evolves the gas for 2 Myr with chemistry on and the driving off.

## What happens

1. **Drive** (`inputs/drive64.athinput`). AthenaK stirs a periodic box,
   32 pc on a side with 64^3 cells, at a mean density n_H = 10 cm^-3. This run
   is hydrodynamics only, starting at T = 10^4 K and running to t = 12 code
   units.
2. **Shared initial condition.** The last snapshot of the driven run is
   written as a legacy VTK file, `shared_ic64.vtk`.
3. **Chemistry runs.** AthenaK (`inputs/turb64_si.athinput`, or
   `turb64_bdf.athinput` for BDF) and Athena++ (`inputs/turb64_cvode.athinput`)
   both read `shared_ic64.vtk` and run for 2 Myr with GOW17 chemistry.
4. **Comparison.** The script plots the history of both runs and slices of the
   final state side by side, and makes a movie.

Both codes use the same hydro scheme (RK2, PLM, HLLC), update the chemistry
once per hydro step, and use a constant H2 formation rate on grains,
3e-17 cm^3 s^-1 (`is_kgrH2_const = true` in AthenaK). The radiation field is
uniform, G0 = 1e-6; there is no six-ray shielding. The `_n1_` inputs repeat the
test at n_H = 1 cm^-3 in an unshielded field, G0 = 1.

## What you need

**AthenaK** from the `semi-implicit` branch of my fork, which has the two
problem generators this test needs (`GOW17_turb` for the driving and
`GOW17_read_vtk` for reading the snapshot):

    git clone --recursive -b semi-implicit https://github.com/jeonggyukim/athenak.git ~/athenak
    cd ~/athenak && git checkout 2deb5e40
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON
    make -C build -j8

**Athena++** at upstream commit `823614c`, with one small patch from this
directory (upstream's `read_vtk.cpp` does not compile with MPI), built with
CVODE from SUNDIALS (version 6.7.0 was used):

    git clone https://github.com/PrincetonUniversity/athena.git ~/athena
    cd ~/athena && git checkout 823614c
    git apply <this repo>/reproduce/patches/athena-pp-read_vtk-mpi.diff
    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> -mpi
    make -j8
    gcc -Wall -O2 -o vis/vtk/join_vtk++ vis/vtk/join_vtk++.c -lm

**Python 3** with numpy and matplotlib, and `ffmpeg` for the movie.

## Running it

From the top of this repository:

    SMOKE=1 bash reproduce/run_turb64.sh    # quick check: every step for 20 cycles, ~2 min
    bash reproduce/run_turb64.sh            # the full test

You can also run single steps by name, for example only the shared initial
condition, or only the chemistry runs and plots once it exists:

    bash reproduce/run_turb64.sh drive ic
    bash reproduce/run_turb64.sh si cvode join plot

The steps and how long each took on an Apple-silicon MacBook Pro with 8 cores:

| step | what it does | time |
| --- | --- | --- |
| `drive` | drive the turbulence in AthenaK (hydro only, 1 rank) | 4.6 min |
| `ic` | write the last snapshot as `shared_ic64.vtk` | seconds |
| `si` | AthenaK chemistry run, semi-implicit | 28 s |
| `bdf` | AthenaK chemistry run, Kokkos BDF (rtol 1e-2, atol 1e-15) | about 35-90 min |
| `cvode` | Athena++ chemistry run, CVODE (rtol 1e-2, atol 1e-15) | 17.5 min |
| `join` | merge Athena++'s per-MeshBlock VTK files | seconds |
| `plot` | histories, slice comparison and movie | about 1 min |

## Where things go

The script looks for the codes in `~/athenak` and `~/athena` and writes all
output to `data/` in this repository, which git ignores. To use other places,
set any of these before running: `ATHENAK_ROOT`, `ATHENAPP_ROOT`, `AKBIN`,
`APBIN`, `CHEMDATA` (output), `PYATHENA`, `NP` (number of MPI ranks). The file
`machines/default.sh` lists them; the other files in `machines/` are my own
machines and are used only there.

Next to its output, each run saves `versions.txt` with the commit of each
code, plus any uncommitted changes, so a result can always be traced back to
the code that made it.

## Notes

- Driving the box again on another machine or compiler gives the same
  `shared_ic64.vtk` (it was identical, bit for bit, with g++-15 serial and
  g++-16 with MPI). Small differences on other CPUs would not matter, because
  both codes read whichever file was made.
- `log.md` records what was run, on which commits, and what was found.
- `inputs/athinput.pdr1d` is a separate Athena++ 1D PDR setup (six-ray, CVODE)
  and is not used by `run_turb64.sh`.
