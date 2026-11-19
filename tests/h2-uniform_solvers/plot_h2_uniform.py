"""H2 uniform test (tst/inputs/H2_uniform_test.athinput): H2, H and internal energy
from forward Euler, Kokkos BDF and the semi-implicit solver against the analytic
solution in tst/test_suite/chemistry/test_H2_uniform_gpu.py.
Runs: ~/Documents/Athenak-chem-tests/h2-uniform-solvers/runs/r_<solver>/ (Mac CPU, branch semi-implicit-base)."""
import sys, glob, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
RUN_DIR = os.path.expanduser('~/Documents/Athenak-chem-tests/h2-uniform-solvers/runs')
src = open(os.path.expanduser(
    '~/Projects/athenak-semi-implicit-base/tst/test_suite/chemistry/test_H2_uniform_gpu.py')).read()
exec(src[src.index('def H2_uniform_analytical_solution'):src.index('def H2_uniform_l1_errors')])
RUNS = [('forward_euler', 'forward Euler', '#009E73'), ('kokkos_BDF', 'Kokkos BDF', '#0072B2'),
        ('semi_implicit', 'semi-implicit', '#D55E00')]
def load(s):
    t, rows = [], []
    for f in sorted(glob.glob(f'{RUN_DIR}/r_{s}/tab/*.tab')):
        L = open(f).read().splitlines(); t.append(float(L[0].split('time=')[1].split()[0]))
        rows.append(dict(zip(L[1].replace('#', '').split(), map(float, L[2].split()))))
    return np.array(t), {k: np.array([r[k] for r in rows]) for k in rows[0]}
fig, ax = plt.subplots(2, 3, figsize=(15, 8), sharex=True)
for s, lab, c in RUNS:
    t, d = load(s)
    u, nH, fH, fH2, _, e = H2_uniform_analytical_solution(t)
    eint = d['eint'] / nH * u
    for j, (sim, ana) in enumerate([(d['s_00_chem_H2'], fH2), (d['s_01_chem_H'], fH), (eint, e)]):
        ax[0, j].plot(t, sim, color=c, lw=1.6, label=lab)
        with np.errstate(all='ignore'):
            ax[1, j].plot(t, sim / ana - 1, color=c, lw=1.4)
t = np.linspace(0, t.max(), 400); u, nH, fH, fH2, _, e = H2_uniform_analytical_solution(t)
for j, (ana, name) in enumerate([(fH2, 'x$_{H_2}$'), (fH, 'x$_H$'), (e, 'internal energy per H [erg]')]):
    ax[0, j].plot(t, ana, 'k--', lw=1.2, label='analytic' if j == 0 else None)
    ax[0, j].set_title(name, fontsize=14); ax[1, j].set_xlabel('t [code units]', fontsize=12)
    ax[1, j].axhline(0, color='k', lw=0.6); ax[1, j].set_ylabel('simulation / analytic - 1', fontsize=11)
ax[0, 2].set_yscale('log')
ax[0, 0].legend(fontsize=12, frameon=False)
for a in ax.ravel(): a.tick_params(labelsize=11); a.grid(alpha=0.25)
fig.suptitle('H2 uniform test: three ODE solvers against the analytic solution (one cell, Mac CPU)', fontsize=14)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = sys.argv[1] if len(sys.argv) > 1 else f'{HERE}/h2_uniform_solvers.png'
fig.savefig(out, dpi=200); print('wrote', out)
