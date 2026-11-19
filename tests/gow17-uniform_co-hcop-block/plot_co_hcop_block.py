"""CO/HCO+ block in the semi-implicit solver: timescales and error with the block on
and off. GOW17 uniform time series (n_H = 109.21 cm^-3, G0 = 1e-6, xi_CR = 2e-16 s^-1).

(a) Timescales along the converged BDF run, from the GOW17 rate coefficients
    (src/chemistry/network/gow17.hpp UpdateRates_, branch semi-implicit 446ad343):
    HCO+ recombination 1/d2, CO -> HCO+ 1/cc, total CO destruction 1/d1, CO net
    change x_CO/|f_CO|, and the hydro step.
(b) Worst CO and HCO+ error over the run against dt_sub_frac (semi_implicit_cfl on
    that branch), block on and off.
(c) x_CO and x_HCO+ against time, reference and dt_sub_frac = 0.1.
(d) Signed error X/X_ref - 1 of (c).

Runs: run.sh and run_h2_limit.sh in this directory.
"""
import argparse
import glob
import os

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SPEC = ["He+", "OHx", "CHx", "CO", "C+", "HCO+", "H2", "H+", "H3+", "H2+", "O+", "Si+"]
IONS = ["He+", "C+", "HCO+", "H3+", "H2+", "H+", "O+", "Si+"]
XC, XO, XHE, GAMMA = 1.6e-4, 3.2e-4, 0.1, 5.0 / 3.0
XI_CR, G0 = 2e-16, 1e-6
LENGTH, MASS, TIME = 3.0856775809623245e18, 6.882615081124725e31, 3.08567758e13
KB, MYR = 1.380649e-16, 3.15576e13
U_E = MASS / (LENGTH * TIME**2)
FRACS = [0.3, 0.1, 0.03, 0.01]


def read_run(d):
    files = sorted(glob.glob(os.path.join(d, "tab", "*.tab")))
    t, rows = [], []
    for f in files:
        with open(f) as fh:
            head = fh.readline()
            cols = fh.readline().split()[1:]
            row = fh.readline().split()
        t.append(float(head.split("time=")[1].split()[0]))
        rows.append([float(v) for v in row])
    a = np.array(rows)
    out = {"t": np.array(t) * TIME / MYR}
    for name in ["dens", "eint"]:
        out[name] = a[:, cols.index(name)]
    for s in SPEC:
        out[s] = a[:, next(i for i, c in enumerate(cols) if c.endswith("chem_" + s))]
    return out


def rates(r):
    n = r["dens"]
    xe = sum(r[s] for s in IONS)
    cv = KB * ((1 - 2 * r["H2"]) + r["H2"] + XHE + xe) / (GAMMA - 1)
    T = np.maximum(r["eint"] * U_E / n / cv, 1.0)
    xC = np.maximum(XC - r["HCO+"] - r["CHx"] - r["CO"] - r["C+"], 0)
    xO = np.maximum(XO - r["HCO+"] - r["OHx"] - r["CO"] - r["O+"], 0)
    d2 = 1.06e-5 * T**-0.64 * n * xe
    cc = 6.52 * XI_CR + 1.7e-9 * n * r["H3+"]
    d1 = (92 * XI_CR * 2 * r["H2"] + 6.52 * XI_CR + 1.7e-9 * n * r["H3+"]
          + 1.6e-9 * n * r["He+"] + 2.4e-10 * G0)
    a = (7.0e-11 * n * r["CHx"] * xO
         + 7.95e-10 * T**-0.339 * np.exp(0.108 / T) * n * r["OHx"] * xC)
    f_co = a + d2 * r["HCO+"] - d1 * r["CO"]
    return dict(T=T, d1=d1, d2=d2, cc=cc, f_co=f_co)


def signed_err(run, ref, s):
    return run[s] / np.interp(run["t"], ref["t"], ref[s]) - 1


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--data-dir", default=".")
    p.add_argument("--out", default="gow17_co_hcop_block.png")
    p.add_argument("--prefix", default="",
                   help="run-directory prefix, e.g. h2lim_ for out_h2lim_f0.1_block_true")
    p.add_argument("--limit-note", default="substep limited by H2, CO and energy")
    args = p.parse_args()
    ref = read_run(os.path.join(args.data_dir, "out_ref"))
    runs = {(f, b): read_run(os.path.join(args.data_dir,
                                          f"out_{args.prefix}f{f}_block_{b}"))
            for f in FRACS for b in ["true", "false"]}
    k = rates(ref)
    t = ref["t"]

    plt.rcParams.update({"font.size": 13})
    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5), constrained_layout=True)
    on, off = "#1f77b4", "#d95f02"

    a = ax[0, 0]
    tm = t > 0
    a.loglog(t[tm], 1 / k["d2"][tm] / MYR, label=r"HCO$^+$ recombination $1/d_2$")
    a.loglog(t[tm], 1 / k["cc"][tm] / MYR, label=r"CO$\to$HCO$^+$ $1/c$")
    a.loglog(t[tm], 1 / k["d1"][tm] / MYR, "--", label=r"CO destruction $1/d_1$")
    a.loglog(t[tm], ref["CO"][tm] / np.abs(k["f_co"][tm]) / MYR,
             label=r"CO net change $x_{\rm CO}/|f_{\rm CO}|$")
    dt = np.diff(t)
    a.loglog(t[1:][dt > 0], dt[dt > 0], "k:", label="hydro step")
    a.set_xlabel("t [Myr]")
    a.set_ylabel("timescale [Myr]")
    a.set_title("(a) timescales, converged BDF run")
    a.legend(fontsize=10, loc="lower right")
    frac = k["cc"][-1] / k["d1"][-1]
    a.text(0.03, 0.95, f"at 30 Myr: T = {k['T'][-1]:.0f} K, "
           f"$c/d_1$ = {frac:.2f}", transform=a.transAxes, va="top", fontsize=11)

    a = ax[0, 1]
    for b, col, lab in [("true", on, "block on"), ("false", off, "block off")]:
        for s, mk, ls in [("CO", "o", "-"), ("HCO+", "s", "--")]:
            e = [np.max(np.abs(signed_err(runs[(f, b)], ref, s))) for f in FRACS]
            a.loglog(FRACS, e, marker=mk, ls=ls, color=col,
                     label=f"{s.replace('+', '$^+$')}, {lab}")
    a.set_xlabel("dt_sub_frac")
    a.set_ylabel(r"max$_t$ $|X/X_{\rm ref} - 1|$")
    a.set_title("(b) worst error over 30 Myr")
    a.legend(fontsize=10)

    a = ax[1, 0]
    for s, ls in [("CO", "-"), ("HCO+", "--")]:
        a.loglog(t[tm], ref[s][tm], "k", ls=ls, lw=2.5, alpha=0.5,
                 label=f"{s.replace('+', '$^+$')}, BDF reference")
        for b, col, lab in [("true", on, "on"), ("false", off, "off")]:
            r = runs[(0.1, b)]
            m = r["t"] > 0
            a.loglog(r["t"][m], r[s][m], color=col, ls=ls, lw=1.2,
                     label=f"{s.replace('+', '$^+$')}, block {lab}")
    a.set_xlabel("t [Myr]")
    a.set_ylabel("abundance per H")
    a.set_title("(c) dt_sub_frac = 0.1")
    a.legend(fontsize=10)

    a = ax[1, 1]
    for s, ls in [("CO", "-"), ("HCO+", "--")]:
        for b, col, lab in [("true", on, "on"), ("false", off, "off")]:
            r = runs[(0.1, b)]
            m = r["t"] > 0
            a.semilogx(r["t"][m], signed_err(r, ref, s)[m], color=col, ls=ls,
                       label=f"{s.replace('+', '$^+$')}, block {lab}")
    a.axhline(0, color="0.5", lw=0.8)
    a.set_xlabel("t [Myr]")
    a.set_ylabel(r"$X/X_{\rm ref} - 1$")
    a.set_title("(d) signed error, dt_sub_frac = 0.1")
    a.legend(fontsize=10)

    fig.suptitle(args.limit_note, fontsize=14)
    fig.savefig(args.out, dpi=150)
    print("wrote", os.path.abspath(args.out))
    for f in FRACS:
        print(f, *[f"{b}:{np.max(np.abs(signed_err(runs[(f, b)], ref, 'CO'))):.3g}"
                   for b in ["true", "false"]])


if __name__ == "__main__":
    main()
