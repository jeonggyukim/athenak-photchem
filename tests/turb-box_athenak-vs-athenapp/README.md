# Decaying turbulence with GOW17 chemistry: AthenaK and Athena++

We compare the chemical and thermal evolution of turbulent gas computed with
the GOW17 network in AthenaK and in Athena++. In AthenaK the rate equations are
integrated with the semi-implicit solver (or with Kokkos BDF); in Athena++
they are integrated with CVODE. Agreement between the two, cell by cell, tests
the AthenaK implementation of the network and of the solver against an
established code.

The two codes generate turbulence with different driving schemes, so a direct
comparison requires a common initial state. We therefore drive turbulence
once, in AthenaK, in a periodic box of side 32 pc resolved by 64^3 cells, at a
mean hydrogen density n_H = 10 cm^-3 and an initial temperature of 10^4 K,
without chemistry. The final snapshot of this run is the initial condition for
both codes, which then follow its evolution for 2 Myr with chemistry included
and the driving switched off. Both codes use the same hydrodynamic scheme (RK2
time integration, PLM reconstruction, HLLC fluxes), update the chemistry once
per hydrodynamic step, and adopt a constant H2 formation rate coefficient on
grains, 3e-17 cm^3 s^-1. The gas is exposed to a uniform, weak radiation field,
G0 = 1e-6 in Draine units, so that its chemistry is set by cosmic rays and
collisions; shielding is not computed. A second set of runs, labeled `_n1_`,
repeats the test at n_H = 1 cm^-3 in an unattenuated field, G0 = 1.

## Reproducing the results

The AthenaK version used here is commit `2deb5e40` on the `semi-implicit`
branch of `github.com/jeonggyukim/athenak`, which contains the problem
generators for driving (`GOW17_turb`) and for reading the shared initial state
(`GOW17_read_vtk`):

    git clone --recursive -b semi-implicit https://github.com/jeonggyukim/athenak.git ~/athenak
    cd ~/athenak && git checkout 2deb5e40
    cmake -B build -DCMAKE_BUILD_TYPE=Release -DAthena_ENABLE_MPI=ON && make -C build -j8

The Athena++ version is upstream commit `823614c`, compiled with CVODE
(SUNDIALS 6.7.0) after a one-line fix that lets `read_vtk.cpp` build with MPI:

    git clone https://github.com/PrincetonUniversity/athena.git ~/athena
    cd ~/athena && git checkout 823614c
    git apply <this repo>/tests/turb-box_athenak-vs-athenapp/patches/athena-pp-read_vtk-mpi.diff
    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> -mpi
    make -j8
    gcc -Wall -O2 -o vis/vtk/join_vtk++ vis/vtk/join_vtk++.c -lm

The figures are made with Python 3 (numpy, matplotlib), and the movie with
`ffmpeg`. From the top of this repository, the whole comparison is rerun with

    bash tests/turb-box_athenak-vs-athenapp/run_turb64.sh

which drives the turbulence (about 5 min on an 8-core laptop), writes the
shared initial state, evolves it with the AthenaK semi-implicit solver (30 s),
Kokkos BDF (35-90 min) and Athena++ CVODE (18 min), and plots the results.
Individual parts can be run by name (`drive`, `ic`, `si`, `bdf`, `cvode`,
`join`, `plot`), and `SMOKE=1` runs every part for 20 cycles as a quick check
of the setup. The resolution is set by `NX` (default 64, i.e. 64^3 cells); for
example, `NX=128` repeats the whole comparison at 128^3. The codes are expected in `~/athenak` and `~/athena` and the
output is written to `data/`; other locations are set through the variables
listed in `machines/default.sh`. The code versions behind each run are saved
with its output, and `log.md` records the runs made so far and what they
showed.
