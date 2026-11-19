#!/usr/bin/env python
"""Compare the GOW17 ordered semi-implicit sweep against the Kokkos BDF reference.

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
    ~/miniconda3/envs/pyathena/bin/python plot_sweep_vs_bdf.py \
        --data-dir /gpfs/jeonggyukim/athenak-chem/runs/acc4 \
        --cost-dir /gpfs/jeonggyukim/athenak/runs/perf64 \
        --gpu a100 --cost-gpu h200/syn104

"""
import argparse, glob, os, re, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
NH, KB, PCGS = 109.21, 1.380649e-16, 2.3426e-14

# accuracy arms: 4^3 (all cells identical -> cell count cannot affect the answer)
ARMS = [('bdf_lg',          'Kokkos BDF (reference)', '#111111', 2.2, '-'),
        ('n3_lg',           'sweep n=3',   '#c6dbef', 1.4, '-'),
        ('n10_h0_m0_b0_lg', 'sweep n=10',  '#6baed6', 1.4, '-'),
        ('n30_lg',          'sweep n=30',  '#2171b5', 1.4, '-'),
        ('n100_lg',         'sweep n=100', '#08306b', 1.4, '-'),
        ('pc10_lg',         'adaptive on H2/CO/$e_{int}$ (cfl=0.1)', '#D55E00', 1.6, '--')]
# cost arms: label, 64^3 log tag, colour. Order is the bar order.
COST_ARMS = [('hydro only (HLLC)', 'hydro', '#777777'),
             ('sweep n=3',         'o3',    '#c6dbef'),
             ('sweep n=10',        'o10',   '#6baed6'),
             ('sweep n=30',        'o30',   '#2171b5'),
             ('sweep n=100',       'o100',  '#08306b'),
             ('Kokkos BDF',        'bdfN',  '#111111')]
ADAPTIVE_MEAN_SUBSTEPS = 1.2

ap = argparse.ArgumentParser()
ap.add_argument('--data-dir', default='.', help='4^3 run tree holding out_<arm>/tab/')
ap.add_argument('--cost-dir', default=None, help='64^3 run tree holding log_<tag>.txt')
ap.add_argument('--gpu', default='', help='GPU the 4^3 arms ran on')
ap.add_argument('--cost-gpu', default='', help='GPU the 64^3 cost arms ran on')
ap.add_argument('--out', default='gow17_sweep_adaptive.png')
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


COST, missing = [], []
if args.cost_dir:
    for label, tag, colour in COST_ARMS:
        c = cost_per_cycle(args.cost_dir, tag)
        if c is None:
            missing.append('%s (log_%s.txt)' % (label, tag))
        else:
            COST.append((label, c, colour))
    fixed = [(int(lab.split('=')[1]), c) for lab, c, _ in COST if lab.startswith('sweep n=')]
    if fixed:
        rate = float(np.mean([c / n for n, c in fixed]))
        COST.insert(1, ('adaptive (full-run mean)',
                        ADAPTIVE_MEAN_SUBSTEPS * rate, '#D55E00'))
for m in missing:
    print('cost panel: no timing for %s' % m)


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

tref,dref = D['bdf_lg']; onref=lambda t,y: np.interp(tref,t,y)
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
a.legend(fontsize=7,frameon=False,loc='lower left')

a=ax[14]
for tag,lab,c,lw,ls in [('pc30_lg','adaptive cfl=0.3','#f4a582',1.3,'-'),
                        ('pc10_lg','adaptive cfl=0.1','#D55E00',1.4,'-')]:
    try: txt=open(os.path.join(args.data_dir,'log_%s.txt'%tag)).read()
    except IOError: continue
    tt=[float(x) for x in re.findall(r'^elapsed=\S+ cycle=\d+ time=(\S+)', txt, re.M)]
    mm=[float(x) for x in re.findall(r'chemistry substeps: total=\S+ mean=(\S+)', txt)]
    k=min(len(tt),len(mm))
    if k: a.plot(tt[:k], mm[:k], color=c, lw=lw, label=lab)
for n,c in [(3,'#c6dbef'),(10,'#6baed6'),(30,'#2171b5'),(100,'#08306b')]:
    a.axhline(n,color=c,lw=1.0,ls=':')
    a.text(25,n*1.08,'fixed n=%d'%n,fontsize=6.5,color=c,ha='right')
a.set_yscale('log'); a.set_ylabel('substeps per hydro step')
a.set_title('adaptive substep count vs time',fontsize=10)
a.legend(fontsize=7,frameon=False,loc='upper right')

a=ax[15]
if COST:
    nm=[x[0] for x in COST]; cv=[x[1] for x in COST]; cc=[x[2] for x in COST]
    a.barh(range(len(nm)),cv,color=cc,height=0.6,edgecolor='#444444',linewidth=0.4)
    a.set_yticks(range(len(nm))); a.set_yticklabels(nm,fontsize=8); a.invert_yaxis()
    a.set_xscale('log'); a.set_xlim(5e-4,400); a.set_xlabel('CPU time per cycle [s]')
    ref = cv[0]
    for i,c in enumerate(cv):
        a.text(c*1.35,i,'%.3g s   (%.0f$\\times$ first bar)'%(c,c/ref),va='center',fontsize=8)
else:
    a.text(0.5,0.5,'no 64$^3$ timings available',ha='center',va='center',fontsize=11)
    a.set_xticks([]); a.set_yticks([])
a.set_title('cost at 64$^3$ = 262,144 cells%s'%(
    ' (%s)'%args.cost_gpu if args.cost_gpu else ''),fontsize=10)

for a in ax[:15]:
    a.set_xlabel('t [Myr]'); a.set_xscale('log'); a.set_xlim(9e-3,30)
    a.grid(alpha=0.25,lw=0.5); a.tick_params(labelsize=8)
ax[15].grid(False); ax[15].tick_params(labelsize=8)
ax[0].legend(fontsize=8,loc='lower left',frameon=False)
prov = 'evolution, error & substep panels: 4$^3$ cells%s' % (' on %s'%args.gpu if args.gpu else '')
if COST:
    prov += ';  cost panel: 64$^3$ cells%s' % (' on %s'%args.cost_gpu if args.cost_gpu else '')
    if missing:
        prov += '  [%d arm(s) omitted: no timing log]' % len(missing)
fig.suptitle(
 'GOW17: ordered semi-implicit sweep vs Kokkos BDF.  n = substeps per hydro step, '
 '$dt_{chem}=dt_{hydro}/n$; fixed-n arms use no per-cell adaptivity\n'
 'reference solver: Kokkos Kernels variable-order BDF (order 1-5, adaptive internal step)\n'
 'atol=10$^{-6}$, rtol=10$^{-3}$, max 10 Newton iterations, min step factor 0.2, '
 'first internal step auto-selected;  error weight $s_i=$atol$+$rtol$|y_i|$, RMS over all 13 equations\n'
 'numerical Jacobian (forward difference, $\\sqrt{\\epsilon}\\max(|y|,\\mathrm{atol})$), '
 'LU with partial pivoting; temperature-only rates hoisted out of the Jacobian\n'
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
fig.tight_layout(rect=[0,0,1,0.935]); fig.savefig(args.out,dpi=300)
print('wrote %s' % args.out)
