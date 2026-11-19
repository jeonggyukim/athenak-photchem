# athenak-photchem

Tests of photochemistry in AthenaK: the chemistry module (GOW17 and H2 networks,
semi-implicit and Kokkos BDF solvers), and later its coupling to radiation and
the ray tracer. Each test has its inputs, run script and plot script here;
simulation output is not committed, the scripts rebuild it.

| directory | holds |
| --- | --- |
| `reproduce/` | 64^3 decaying turbulence, AthenaK against Athena++ CVODE from a shared initial condition; see [`reproduce/README.md`](reproduce/README.md) |
| `tests/<problem>/` | one folder per test problem: inputs, `run.sh`, plot script, figures in `fig_<name>/` |
| `scripts/` | plotting and comparison scripts shared by the turbulence test and the GOW17 network figures |
| `figures/`, `movie/` | figures and movies of the turbulence test and the GOW17 network |

## Tests

| test | compares |
| --- | --- |
| [`reproduce/`](reproduce/) | turbulent box: AthenaK semi-implicit and BDF against Athena++ CVODE, cell by cell |
| [`tests/gow17-uniform_check-times/`](tests/gow17-uniform_check-times/) | GOW17 uniform box, four conditions: semi-implicit and default BDF against converged BDF at fixed times (prototype of the AthenaK regression test) |
| [`tests/gow17-uniform_co-hcop-block/`](tests/gow17-uniform_co-hcop-block/) | CO/HCO+ solved together or separately, timescales and errors |
| [`tests/gow17-uniform_refresh-rates/`](tests/gow17-uniform_refresh-rates/) | semi-implicit `refresh_rates` and `nsub_max`: accuracy and cost |
| [`tests/gow17-uniform_semi-implicit-vs-bdf/`](tests/gow17-uniform_semi-implicit-vs-bdf/) | semi-implicit with and without error control against converged BDF |
| [`tests/gow17-uniform_bdf-atol/`](tests/gow17-uniform_bdf-atol/) | Kokkos BDF absolute tolerance sweep, accuracy and GPU cost |
| [`tests/gow17-uniform_semi-implicit-species-vs-bdf/`](tests/gow17-uniform_semi-implicit-species-vs-bdf/) | early semi-implicit species evolution against BDF |
| [`tests/h2-uniform_solvers/`](tests/h2-uniform_solvers/), [`tests/h2-advection_solvers/`](tests/h2-advection_solvers/) | H2 network, ODE solvers against the analytic solution |

Each `run.sh` writes its output outside the repository (by default
`~/Documents/Athenak-chem-tests/<study>/`; set `OUT` to change it) and states
the AthenaK commit it was run with.
