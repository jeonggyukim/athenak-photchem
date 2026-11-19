#!/usr/bin/env python
"""GOW17 uniform time series: Kokkos BDF at rtol 1e-3 with atol 1e-6 ... 1e-14
and the semi-implicit solver, against converged BDF (rtol 1e-7, atol 1e-16).

Same layout as tests/gow17-uniform_semi-implicit-species-vs-bdf/plot_sweep_vs_bdf.py.
Runs from job.sbatch in this folder (syntax job on one H200). The last two panels
are the 64^3 cost in two windows: cycles 0-40 and 40 cycles from t = 10 Myr.

    python plot_bdf_atol.py --data-dir runs/acc --cost-dir runs/costA \
        --cost-dir-b runs/costB --gpu H200 --cost-gpu H200

Reads AthenaK tab dumps from `<data-dir>/out_<arm>/tab/` and produces a 4x4 panel
figure: temperature, twelve species, worst-species error vs the BDF, adaptive
substep count vs time, and cost.

Two meshes, deliberately. Accuracy is measured at 4^3 = 64 cells: the pgen is a
uniform static box with vx = 0, so every cell is identical and cell count cannot
change the answer -- the small mesh only makes a full BDF reference to 30 Myr
affordable (about a minute rather than hours). Cost is measured separately at
64^3, because 4^3 timings are latency-bound and say nothing about throughput.

Accuracy is hardware-independent, so the 4^3 arms may be run on any GPU. Cost is
not, so --cost-gpu is labelled separately from --gpu and the two need not match.

Cost per cycle is read from the 64^3 run logs as `cpu time used` divided by
`MeshBlock-cycles`. An arm whose log is missing, or whose run did not reach the
timing line, is dropped from the panel and named on stdout. It is never replaced
by a literal: an earlier version of this script carried hard-coded values for the
hydro-only and BDF bars that no log in the tree reproduced.

The adaptive arm's cost is the one reconstruction. Its 64^3 window is 20 cycles
of pure transient, where the controller is supposed to burst, so timing it there
reports ~n=10 instead of the full-run mean. Substep count is a per-cell physics
property and therefore identical at both meshes, so the full-run cost is taken as
(mean substeps) x (cost per substep), the latter derived from the fixed-n arms.

Run with the pyathena interpreter:
    python plot_semi_implicit_vs_bdf.py \
        --data-dir ~/Documents/Athenak-chem-tests/gow17-semi-implicit-vs-bdf \
        --gpu "Mac CPU" --out gow17_semi_implicit_vs_bdf.png

"""
import argparse, glob, os, re, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
NH, KB, PCGS = 109.21, 1.380649e-16, 2.3426e-14

# accuracy arms: 4^3 (all cells identical -> cell count cannot affect the answer)
ARMS = [('bdf_ref', 'BDF converged (reference)', '#111111', 2.2, '-'),
        ('bdf_a6',  'BDF atol 1e-6 (default)', '#bdbdbd', 1.6, '--'),
        ('bdf_a8',  'BDF atol 1e-8',  '#c6dbef', 1.4, '-'),
        ('bdf_a10', 'BDF atol 1e-10', '#6baed6', 1.4, '-'),
        ('bdf_a12', 'BDF atol 1e-12', '#2171b5', 1.4, '-'),
        ('bdf_a14', 'BDF atol 1e-14', '#08306b', 1.4, '-'),
        ('si',      'semi-implicit, cfl=0.1', '#D55E00', 0.8, '-')]
# cost arms: label, 64^3 log tag, colour. Order is the bar order.
COST_ARMS = [('semi-implicit',   'si',      '#D55E00'),
             ('BDF atol 1e-6',   'bdf_a6',  '#bdbdbd'),
             ('BDF atol 1e-8',   'bdf_a8',  '#c6dbef'),
             ('BDF atol 1e-10',  'bdf_a10', '#6baed6'),
             ('BDF atol 1e-12',  'bdf_a12', '#2171b5'),
             ('BDF atol 1e-14',  'bdf_a14', '#08306b')]
ap = argparse.ArgumentParser()
ap.add_argument('--data-dir', default='.', help='4^3 run tree holding out_<arm>/tab/')
ap.add_argument('--cost-dir', default=None, help='64^3 window A (cycles 0-40): log_<tag>.txt')
ap.add_argument('--cost-dir-b', default=None, help='64^3 window B (40 cycles from 10 Myr)')
ap.add_argument('--gpu', default='', help='GPU the 4^3 arms ran on')
ap.add_argument('--cost-gpu', default='', help='GPU the 64^3 cost arms ran on')
ap.add_argument('--out', default='gow17_bdf_atol.png')
args = ap.parse_args()


def cost_per_cycle(cost_dir, tag):
    """cpu seconds per cycle from an AthenaK log, or None when unavailable."""
    try:
        txt = open(os.path.join(cost_dir, 'log_%s.txt' % tag)).read()
    except IOError:
        return None
    m_cpu = re.search(r'cpu time used\s*=\s*(\S+)', txt)
    m_cyc = re.search(r'MeshBlock-cycles\s*=\s*(\d+)', txt)
    if not (m_cpu and m_cyc and int(m_cyc.group(1)) > 0):
        return None
    return float(m_cpu.group(1)) / int(m_cyc.group(1))


def collect(cost_dir):
    cost, missing = [], []
    if cost_dir:
        for label, tag, colour in COST_ARMS:
            c = cost_per_cycle(cost_dir, tag)
            if c is None:
                missing.append('%s (log_%s.txt)' % (label, tag))
            else:
                cost.append((label, c, colour))
    for m in missing:
        print('cost panel %s: no timing for %s' % (cost_dir, m))
    return cost, missing
COST, missing = collect(args.cost_dir)
COST_B, missing_b = collect(args.cost_dir_b)


def load(tag):
    t, rows = [], []
    pat = os.path.join(args.data_dir, 'out_%s/tab/%s.hydro_w.*.tab' % (tag, tag))
    for f in sorted(glob.glob(pat)):
        L = open(f).read().splitlines()
        t.append(float(L[0].split('time=')[1].split()[0]))
        rows.append(dict(zip(L[1].replace('#','').split(), [float(v) for v in L[2].split()])))
    if not rows:
        raise SystemExit('no tab files matched %s' % pat)
    return np.array(t), {k: np.array([r[k] for r in rows]) for k in rows[0]}
D = {a[0]: load(a[0]) for a in ARMS}
SPEC = [('s_06_chem_H2','H$_2$'),('s_07_chem_H+','H$^+$'),('s_00_chem_He+','He$^+$'),
        ('s_04_chem_C+','C$^+$'),('s_03_chem_CO','CO'),('s_02_chem_CHx','CH$_x$'),
        ('s_01_chem_OHx','OH$_x$'),('s_05_chem_HCO+','HCO$^+$'),('s_08_chem_H3+','H$_3^+$'),
        ('s_09_chem_H2+','H$_2^+$'),('s_10_chem_O+','O$^+$'),('s_11_chem_Si+','Si$^+$')]
def temp(d):
    xe = sum(d['s_%02d_chem_%s'%(i,n)] for i,n in
             [(0,'He+'),(4,'C+'),(5,'HCO+'),(8,'H3+'),(9,'H2+'),(7,'H+'),(10,'O+'),(11,'Si+')])
    return d['eint']*PCGS/NH/(1.5*KB*(1.1-d['s_06_chem_H2']+xe))

fig, ax = plt.subplots(4,4,figsize=(17,15.5)); ax=ax.ravel()
for tag,lab,c,lw,ls in ARMS:
    t,d = D[tag]; ax[0].plot(t, temp(d), color=c, lw=lw, ls=ls, label=lab)
ax[0].set_ylabel('T [K]'); ax[0].set_yscale('log'); ax[0].set_title('temperature', fontsize=10)
for i,(k,lab) in enumerate(SPEC):
    a=ax[i+1]
    for tag,_,c,lw,ls in ARMS:
        t,d=D[tag]; a.plot(t,d[k],color=c,lw=lw,ls=ls)
    a.set_yscale('log'); a.set_title(lab,fontsize=10); a.set_ylabel('$x_i$ (per H)',fontsize=8)

tref,dref = D['bdf_ref']; onref=lambda t,y: np.interp(tref,t,y)
a=ax[13]
for tag,lab,c,lw,ls in ARMS[1:]:
    t,d=D[tag]; e=np.zeros_like(tref)
    for k,_ in SPEC+[('eint','eint')]:
        with np.errstate(divide='ignore',invalid='ignore'):
            e=np.maximum(e,np.nan_to_num(np.abs(onref(t,d[k])/dref[k]-1.0)))
    a.plot(tref,e,color=c,lw=lw,label=lab)
    # temperature is plotted separately: a loose rtol degrades T, which feeds
    # every rate coefficient, and folding it into a max over 13 equations hides
    # which one moved.
    with np.errstate(divide='ignore',invalid='ignore'):
        eT=np.nan_to_num(np.abs(onref(t,temp(d))/temp(dref)-1.0))
    a.plot(tref,eT,color=c,lw=0.9,ls=':')
a.set_yscale('log'); a.set_ylabel('max$_i$ |$x_i/x_i^{BDF}-1$|')
a.set_title('worst-species error vs BDF (dotted: T)', fontsize=10)
a.axhline(0.01,color='#999999',lw=0.8); a.text(0.011,0.012,'1%',fontsize=8,color='#666666')

def draw_cost(a, cost, window):
    if cost:
        nm=[x[0] for x in cost]; cv=[x[1] for x in cost]; cc=[x[2] for x in cost]
        a.barh(range(len(nm)),cv,color=cc,height=0.6,edgecolor='#444444',linewidth=0.4)
        a.set_yticks(range(len(nm))); a.set_yticklabels(nm,fontsize=8); a.invert_yaxis()
        a.set_xscale('log'); a.set_xlim(min(cv)/2, max(cv)*150)
        a.set_xlabel('wall time per cycle [s]')
        ref = cv[0]
        for i,c in enumerate(cv):
            a.text(c*1.35,i,'%.3g s   (%.1f$\\times$)'%(c,c/ref),va='center',fontsize=8)
    else:
        a.text(0.5,0.5,'no 64$^3$ timings available',ha='center',va='center',fontsize=11)
        a.set_xticks([]); a.set_yticks([])
    a.set_title('64$^3$ cost, %s' % window, fontsize=10)
    a.grid(False); a.tick_params(labelsize=8)
draw_cost(ax[14], COST, 'cycles 0-40 (t < 1 Myr)')
draw_cost(ax[15], COST_B, '40 cycles from 10 Myr')

for a in ax[:14]:
    a.set_xlabel('t [Myr]'); a.set_xscale('log'); a.set_xlim(9e-3,30)
    a.grid(alpha=0.25,lw=0.5); a.tick_params(labelsize=8)
handles, labels = ax[0].get_legend_handles_labels()
prov = 'evolution & error panels: 4$^3$ cells%s' % (' on %s'%args.gpu if args.gpu else '')
if COST or COST_B:
    prov += ';  cost panels: 64$^3$ cells, box scaled so $\\Delta x$ and $dt$ match%s' % (
        ' on %s'%args.cost_gpu if args.cost_gpu else '')
    if missing or missing_b:
        prov += '  [%d arm(s) omitted: no timing log]' % (len(missing) + len(missing_b))
fig.suptitle(
 'GOW17: Kokkos BDF absolute tolerance.  Reference: Kokkos BDF converged, rtol=10$^{-7}$, atol=10$^{-16}$\n'
 'BDF runs: rtol=10$^{-3}$ with atol from 10$^{-6}$ (the default) to 10$^{-14}$; '
 'semi-implicit: semi_implicit_cfl=0.1, nsub_max=100 (branch semi-implicit-base 9731520a)\n'
 'uniform static box (no flow), $n_H$=109.21 cm$^{-3}$, $\\mu_H$=1.4, $T(t{=}0)$=154 K, '
 'all species $x_i(t{=}0)=10^{-6}$\n'
 '$Z_d=Z_g=1$, $x_{He}$=0.1, $x_C$=1.6$\\times$10$^{-4}$, '
 '$x_O$=3.2$\\times$10$^{-4}$, $x_{Si}$=1.7$\\times$10$^{-6}$;  '
 'radiation: $G_0=10^{-6}$ Draine (fully shielded, no FUV), $\\zeta_{CR}=2\\times10^{-16}$ s$^{-1}$/H\n'
 '$\\Delta x$=0.0156 pc, box 0.0625 pc; $dt_{hydro}$ set by the hydro CFL: '
 '9.68$\\times$10$^{-3}$ to 3.52$\\times$10$^{-2}$ Myr, median 2.83$\\times$10$^{-2}$.  '
 'CFL 0.8, rk2, $\\gamma$=5/3, $t_{lim}$=30 Myr (1033 cycles)\n'
 + prov,
 fontsize=13, y=0.995, va='top')
for a in ax:
    a.title.set_fontsize(14)
    a.xaxis.label.set_fontsize(12); a.yaxis.label.set_fontsize(12)
    a.tick_params(labelsize=11)
for a in ax[14:]:
    for txt in a.texts: txt.set_fontsize(11)
fig.subplots_adjust(left=0.05, right=0.985, bottom=0.04, top=0.83, hspace=0.42, wspace=0.42)
leg = fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.905), ncol=4,
                 frameon=False, fontsize=14, handlelength=3.5, columnspacing=3.0)
for h in leg.legend_handles:
    h.set_linewidth(3.0)
fig.savefig(args.out,dpi=300)
print('wrote %s' % args.out)
