# Decaying turbulence with GOW17 chemistry: AthenaK and Athena++

This test evolves the same turbulent gas with the GOW17 chemical network in
AthenaK and in Athena++ and compares the two solutions cell by cell. AthenaK
integrates the chemistry with its semi-implicit solver or with Kokkos BDF;
Athena++ uses CVODE.

Because the two codes drive turbulence differently, the turbulence is generated
once, in AthenaK, and both codes start from the same snapshot. AthenaK first
drives a periodic box 32 pc on a side, resolved by 64^3 cells, at a mean
hydrogen density n_H = 10 cm^-3, without chemistry. Its final state is written
as `shared_ic64.vtk`. Each code then reads this file and evolves it for 2 Myr
with chemistry and without driving. The radiation field is uniform, with
G0 = 1e-6 in Draine units, and no shielding is computed. The inputs labeled
`_n1_` repeat the test at n_H = 1 cm^-3 in an unattenuated field, G0 = 1.

## Codes

AthenaK is taken from the `semi-implicit` branch of the fork below, which
contains the problem generators `GOW17_turb` (driving) and `GOW17_read_vtk`
(reading the snapshot):

    git clone --recursive -b semi-implicit https://github.com/jeonggyukim/athenak.git ~/athenak
    cd ~/athenak && git checkout 2deb5e40
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON && make -C build -j8

Athena++ is taken at commit `823614c`, with the patch in `patches/` that lets
`read_vtk.cpp` compile with MPI, and with CVODE from SUNDIALS 6.7.0:

    git clone https://github.com/PrincetonUniversity/athena.git ~/athena
    cd ~/athena && git checkout 823614c
    git apply <this repo>/reproduce/patches/athena-pp-read_vtk-mpi.diff
    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> -mpi
    make -j8
    gcc -Wall -O2 -o vis/vtk/join_vtk++ vis/vtk/join_vtk++.c -lm

The analysis requires Python 3 with numpy and matplotlib, and `ffmpeg` for the
movie.

## Running the test

The script `run_turb64.sh` performs each step in turn and is run from the top
of this repository. A short run of every step (20 cycles each, about 2 min)
checks the setup before the full test:

    SMOKE=1 bash reproduce/run_turb64.sh    # check the setup
    bash reproduce/run_turb64.sh            # full test
    bash reproduce/run_turb64.sh si plot    # selected steps only

| step | action | time on an 8-core laptop |
| --- | --- | --- |
| `drive` | drive the turbulence in AthenaK | 4.6 min |
| `ic` | write `shared_ic64.vtk` | seconds |
| `si`, `bdf` | AthenaK chemistry, semi-implicit or Kokkos BDF | 28 s; 35-90 min |
| `cvode` | Athena++ chemistry with CVODE | 17.5 min |
| `join`, `plot` | merge the Athena++ output; plot histories, slices and a movie | about 1 min |

By default the script looks for the codes in `~/athenak` and `~/athena` and
writes its output to `data/`, which git ignores. The variables `ATHENAK_ROOT`,
`ATHENAPP_ROOT`, `CHEMDATA` (output directory), `NP` (number of MPI ranks) and
`PYTHON` (interpreter, `python3` by default) override these choices; all are
listed in `machines/default.sh`. Each run records the code versions it used in
`versions.txt`, and `log.md` lists earlier runs and their results.
