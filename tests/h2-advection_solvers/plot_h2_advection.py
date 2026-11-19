"""H2 advection test (tst/inputs/H2_advection_test.athinput): a Gaussian H profile
advected one box crossing to t = 5 pc/(km/s) at nx1 = 32-256, for forward Euler,
Kokkos BDF and the semi-implicit solver, against the analytic solution in
tst/test_suite/chemistry/test_H2_advection_gpu.py.
Left: final x_H profiles at nx1 = 64. Middle: their difference from the analytic
profile. Right: L1 error of x_H2 against nx1, with the test's limits.
Runs: ~/Documents/Athenak-chem-tests/h2-advection-solvers/runs/r_<solver>_<nx1>/ (Mac CPU, branch semi-implicit-base)."""
import sys, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker
HERE = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.path.expanduser('~/Documents/Athenak-chem-tests/h2-advection-solvers/runs')
REPO = os.path.expanduser('~/Projects/athenak-semi-implicit-base')
sys.path.insert(0, REPO + '/vis/python'); import athena_read
src = open(REPO + '/tst/test_suite/chemistry/test_H2_advection_gpu.py').read()
exec(src[src.index('class constants'):src.index('def assert_constant')])
T = 5 * constants.pc_cgs / constants.km_s_cgs
NX = [32, 64, 128, 256]
LIM = 1.1 * np.array([0.008902357671992878, 0.0033214914276244742,
                      0.001289147168499853, 0.0007298300086431499])
RUNS = [('forward_euler', 'forward Euler', '#009E73', 's'), ('kokkos_BDF', 'Kokkos BDF', '#0072B2', 'o'),
        ('semi_implicit', 'semi-implicit', '#D55E00', '^')]
def final(s, n):
    return athena_read.tab(f'{RUN_DIR}/r_{s}_{n}/tab/H2_advection.hydro_w.00001.tab')
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
for s, lab, c, m in RUNS:
    d = final(s, 64); x = d['x1v']; fH, fH2 = H2_advection_analytical_solutions(x, T, 1.5, 0.1)
    ax[0].plot(x, d['s_01_chem_H'], color=c, lw=1.6, label=lab)
    ax[1].plot(x, d['s_01_chem_H'] - fH, color=c, lw=1.6)
    l1 = []
    for n in NX:
        d = final(s, n); fH, fH2 = H2_advection_analytical_solutions(d['x1v'], T, 1.5, 0.1)
        l1.append(np.mean(np.abs(fH2 - d['s_00_chem_H2'])))
    ax[2].loglog(NX, l1, color=c, marker=m, lw=1.4, ms=7, label=lab)
xf = np.linspace(x.min(), x.max(), 400); fH, _ = H2_advection_analytical_solutions(xf, T, 1.5, 0.1)
ax[0].plot(xf, fH, 'k--', lw=1.2, label='analytic')
ax[2].loglog(NX, LIM, 'k:', lw=1.4, label='test limit')
ax[0].set_xlabel('x [pc]', fontsize=12); ax[0].set_ylabel('x$_H$ at t = 5 pc/(km/s)', fontsize=12)
ax[0].set_title('x$_H$ profile, nx1 = 64', fontsize=14)
ax[1].set_xlabel('x [pc]', fontsize=12); ax[1].set_ylabel('x$_H$ - analytic', fontsize=12)
ax[1].set_title('difference from analytic, nx1 = 64', fontsize=14); ax[1].axhline(0, color='k', lw=0.6)
ax[2].set_xlabel('nx1', fontsize=12); ax[2].set_ylabel('L1 error of x$_{H_2}$', fontsize=12)
ax[2].set_title('convergence', fontsize=14); ax[2].set_xticks(NX); ax[2].set_xticklabels(NX); ax[2].xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
ax[0].legend(fontsize=11, frameon=False); ax[2].legend(fontsize=11, frameon=False)
for a in ax: a.tick_params(labelsize=11); a.grid(alpha=0.25)
fig.suptitle('H2 advection test: three ODE solvers against the analytic solution (Mac CPU)', fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.93))
out = sys.argv[1] if len(sys.argv) > 1 else f'{HERE}/h2_advection_solvers.png'
fig.savefig(out, dpi=200); print('wrote', out)
