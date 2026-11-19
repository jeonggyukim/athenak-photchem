# Decaying turbulence with GOW17 chemistry: AthenaK vs Athena++

The same 64^3 turbulent box with GOW17 chemistry, run in AthenaK
(semi-implicit or Kokkos BDF) and in Athena++ (CVODE), compared cell by cell.

AthenaK drives the turbulence (32 pc box, n_H = 10 cm^-3, hydro only), and its
last snapshot is saved as `shared_ic64.vtk`. Both codes start from that file and
evolve it for 2 Myr with chemistry on and the driving off. Radiation is a
uniform G0 = 1e-6 (no six-ray shielding); the `_n1_` inputs repeat the test at
n_H = 1 cm^-3 and G0 = 1.

## Codes

AthenaK, `semi-implicit` branch of my fork (it has the `GOW17_turb` and
`GOW17_read_vtk` problem generators):

    git clone --recursive -b semi-implicit https://github.com/jeonggyukim/athenak.git ~/athenak
    cd ~/athenak && git checkout 2deb5e40
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON && make -C build -j8

Athena++ at `823614c`, with the MPI fix in `patches/` and CVODE (SUNDIALS 6.7.0):

    git clone https://github.com/PrincetonUniversity/athena.git ~/athena
    cd ~/athena && git checkout 823614c
    git apply <this repo>/reproduce/patches/athena-pp-read_vtk-mpi.diff
    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> -mpi
    make -j8
    gcc -Wall -O2 -o vis/vtk/join_vtk++ vis/vtk/join_vtk++.c -lm

Python 3 with numpy and matplotlib, and `ffmpeg` for the movie.

## Run

From the top of this repository:

    SMOKE=1 bash reproduce/run_turb64.sh    # every step for 20 cycles, ~2 min
    bash reproduce/run_turb64.sh            # full test
    bash reproduce/run_turb64.sh si plot    # only some steps

| step | does | time (8-core laptop) |
| --- | --- | --- |
| `drive` | drive the turbulence in AthenaK | 4.6 min |
| `ic` | write `shared_ic64.vtk` | seconds |
| `si` / `bdf` | AthenaK chemistry, semi-implicit / Kokkos BDF | 28 s / 35-90 min |
| `cvode` | Athena++ chemistry, CVODE | 17.5 min |
| `join`, `plot` | merge Athena++ output; histories, slices, movie | about 1 min |

The script expects the codes in `~/athenak` and `~/athena` and writes to
`data/` (ignored by git). Override with `ATHENAK_ROOT`, `ATHENAPP_ROOT`,
`CHEMDATA` (output), `NP` (MPI ranks), or `PYATHENA` (the Python interpreter,
`python3` by default); see `machines/default.sh`. Each run saves the code
versions it used in `versions.txt`. `log.md` records past runs and results.
