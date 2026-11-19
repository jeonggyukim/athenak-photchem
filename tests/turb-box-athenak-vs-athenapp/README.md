# Decaying turbulence with GOW17 chemistry: AthenaK and Athena++

We compare the chemical and thermal evolution of turbulent gas computed with
the GOW17 network in AthenaK and in Athena++. In AthenaK the rate equations are
integrated with the semi-implicit solver (or with Kokkos BDF); in Athena++
they are integrated with CVODE. Agreement between the two, cell by cell, tests
the AthenaK implementation of the network and of the solver against an
established code.

The two codes generate turbulence with different driving schemes, so a direct
comparison requires a common initial state. We therefore drive turbulence
once, in AthenaK, in a periodic box of side 32 pc (64^3 cells by default), at a
mean hydrogen density n_H = 10 cm^-3 and an initial temperature of 10^4 K,
without chemistry. The final snapshot of this run is the initial condition for
both codes, which then follow its evolution for 2 Myr with chemistry included
and the driving switched off. Both codes use the same hydrodynamic scheme (RK2
time integration, PLM reconstruction, HLLC fluxes), update the chemistry once
per hydrodynamic step, and adopt a constant H2 formation rate coefficient on
grains, 3e-17 cm^3 s^-1. The gas is exposed to a uniform, weak radiation field,
G0 = 1e-6 in Draine units; shielding is not computed. A second case repeats the
test with the same snapshot scaled to n_H = 1 cm^-3, in an unattenuated field,
G0 = 1. Kokkos BDF and CVODE use the same tolerances (rtol = 1e-2,
atol = 1e-15), so that their runs differ only in the code.

## Reproducing the results

AthenaK must be built with MPI and include the GOW17 chemistry module, the
semi-implicit solver, and the problem generators for driving (`GOW17_turb`) and
for reading the shared initial state (`GOW17_read_vtk`). Until these problem
generators reach the main chemistry branch, they are available in
`github.com/jeonggyukim/athenak` (branch `semi-implicit`).

Athena++ must be configured with the `read_vtk` problem generator, GOW17
chemistry, constant radiation, CVODE and MPI, with HLLC fluxes:

    python configure.py --prob=read_vtk --chemistry=gow17 --chem_radiation=const \
        --chem_ode_solver=cvode --cvode_path=<sundials prefix> --flux=hllc -mpi

Upstream `read_vtk.cpp` does not compile with MPI (it assigns to an undeclared
variable); `patches/athena-pp-read_vtk-mpi.diff` fixes this. The utility
`vis/vtk/join_vtk++` is also needed to merge the Athena++ output.

The figures are made with Python 3 (numpy, matplotlib), and the movie with
`ffmpeg`. From the top of this repository, the whole comparison is rerun with

    bash tests/turb-box-athenak-vs-athenapp/run.sh

which drives the turbulence (about 5 min on an 8-core laptop at 64^3), writes
the shared initial state, evolves it with the AthenaK semi-implicit solver
(30 s), Kokkos BDF (35-90 min) and Athena++ CVODE (18 min), and plots the
results. Individual parts can be run by name (`drive`, `ic`, `si`, `bdf`,
`cvode`, `join`, `plot`), and `SMOKE=1` runs every part for 20 cycles as a quick
check of the setup. `NX` sets the resolution (default 64) and `CASE=n1` selects
the second case. The codes are expected in `~/athenak` and `~/athena` and the
output is written to `data/`; other locations are set through the variables
listed in `machines/default.sh`. The code versions behind each run are saved
with its output, and `log.md` records the runs made so far (AthenaK
`2deb5e40`, Athena++ `823614c`) and what they showed.
