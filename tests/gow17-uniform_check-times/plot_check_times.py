"""Check-time test of the GOW17 uniform problem, one row per condition.

(a) H2, CO, H+ (species that set the substep, and H+); (b) the remaining species;
(c) temperature: converged Kokkos BDF (solid) and semi_implicit at its defaults
(dashed), dense output; vertical lines mark the check times.
(d) error at each check time from the run.sh check runs (check_errors.py
definition): worst species and internal energy, default BDF and semi_implicit.

Dense runs: current case in ~/Documents/Athenak-chem-tests/gow17-uniform_check-times/
dense/, conditions B-D in .../gow17-conditions-sweep/out/ (same binary f133c08e).
Usage: python plot_check_times.py [--out FILE]
"""
import argparse
import glob
import os
import sys

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from check_errors import errors, final  # noqa: E402

DOC = os.path.expanduser("~/Documents/Athenak-chem-tests")
CHECK = f"{DOC}/gow17-uniform_check-times/out"
SW = f"{DOC}/gow17-conditions-sweep/out"
CASES = [
    ("cur", "all species 1e-6, n_H = 109, shielded", ["0.1", "0.4", "3", "30"],
     f"{DOC}/gow17-uniform_check-times/dense/cur_ref",
     f"{DOC}/gow17-uniform_check-times/dense/cur_si"),
    ("B", "molecular start, n_H = 1e2, G0 = 1", ["1e-4", "1e-3", "0.01", "0.1"],
     f"{SW}/molecular_n1e2_cr2e-16_G1_new_ref", f"{SW}/molecular_n1e2_cr2e-16_G1_new100k_si"),
    ("C", "atomic start, n_H = 1e4, shielded", ["0.02", "0.06", "0.15", "1"],
     f"{SW}/atomic_n1e4_cr2e-16_G1e-6_new_ref", f"{SW}/atomic_n1e4_cr2e-16_G1e-6_new100k_si"),
    ("D", r"atomic start, n_H = 1e3, $\xi_{\rm CR}$ = 2e-15", ["0.5", "0.8", "1.1", "5"],
     f"{SW}/atomic_n1e3_cr2e-15_G1e-6_new_ref", f"{SW}/atomic_n1e3_cr2e-15_G1e-6_new100k_si"),
]
MAIN = [("H2", "#1f77b4"), ("CO", "#2ca02c"), ("H+", "#bc8f00")]
REST = [("He+", "#1f77b4"), ("C+", "#ff7f0e"), ("Si+", "#2ca02c"), ("OHx", "#d62728"),
        ("CHx", "#9467bd"), ("HCO+", "#8c564b"), ("H3+", "#e377c2"), ("H2+", "#7f7f7f"),
        ("O+", "#bcbd22")]
IONS = ["He+", "C+", "HCO+", "H3+", "H2+", "H+", "O+", "Si+"]
LENGTH, MASS, TIME = 3.0856775809623245e18, 6.882615081124725e31, 3.08567758e13
KB, MYR, XHE, GAMMA = 1.380649e-16, 3.15576e13, 0.1, 5.0 / 3.0
U_E = MASS / (LENGTH * TIME**2)


def series(d):
    t, rows, cols = [], [], []
    for f in sorted(glob.glob(os.path.join(d, "tab", "*.tab"))):
        with open(f) as fh:
            t.append(float(fh.readline().split("time=")[1].split()[0]))
            cols = fh.readline().split()[1:]
            rows.append([float(v) for v in fh.readline().split()])
    a = np.array(rows)
    x = {c.split("chem_")[-1]: a[:, i] for i, c in enumerate(cols)}
    xe = sum(x[s] for s in IONS)
    cv = KB * ((1 - 2 * x["H2"]) + x["H2"] + XHE + xe) / (GAMMA - 1)
    x["T"] = x["eint"] * U_E / x["dens"] / cv
    return np.array(t) * TIME / MYR, x


def with_transient(t, x, kind):
    """Prepend the condition-B transient (run_B_transient.sh end states, inside the
    first hydro step) to a dense series."""
    ends = sorted(glob.glob(f"{DOC}/gow17-uniform_check-times/transient_B/t*_{kind}"),
                  key=lambda d: float(os.path.basename(d)[1:].rsplit("_", 1)[0]))
    tt, xx = zip(*[series(d) for d in ends])
    tt = np.array([ti[-1] for ti in tt])
    m = t > tt[-1]
    return (np.concatenate([tt, t[m]]),
            {k: np.concatenate([np.array([xi[k][-1] for xi in xx]), x[k][m]]) for k in x})


TEX = {"H2": r"H$_2$", "CO": "CO", "H+": r"H$^+$", "He+": r"He$^+$", "C+": r"C$^+$",
       "Si+": r"Si$^+$", "OHx": r"OH$_x$", "CHx": r"CH$_x$", "HCO+": r"HCO$^+$",
       "H3+": r"H$_3^+$", "H2+": r"H$_2^+$", "O+": r"O$^+$"}


def tex(s):
    return TEX[s]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", default=os.path.join(HERE, "fig_check-times",
                                                  "gow17_check_times.png"))
    args = p.parse_args()
    plt.rcParams.update({"font.size": 16})
    fig, ax = plt.subplots(len(CASES), 4, figsize=(22, 5.0 * len(CASES)),
                           constrained_layout=True)
    for i, (name, title, times, dref, dsi) in enumerate(CASES):
        (tr, r), (ts, s) = series(dref), series(dsi)
        if name == "B":
            tr, r = with_transient(tr, r, "ref")
            ts, s = with_transient(ts, s, "si")
        mr, ms = tr > 0, ts > 0
        for j, group in enumerate([MAIN, REST]):
            a = ax[i, j]
            for k, c in group:
                a.loglog(tr[mr], r[k][mr], color=c, lw=2, label=tex(k))
                a.loglog(ts[ms], s[k][ms], color=c, lw=1.2, ls="--")
            a.set_ylim(bottom=1e-15)
            a.legend(fontsize=11, ncol=3 if j else 1, loc="lower right")
        a = ax[i, 2]
        a.loglog(tr[mr], r["T"][mr], "k", lw=2, label="converged BDF")
        a.loglog(ts[ms], s["T"][ms], color="#D55E00", lw=1.2, ls="--", label="semi_implicit")
        a.legend(fontsize=11)
        for j in range(3):
            for tc in times:
                ax[i, j].axvline(float(tc), color="0.75", lw=0.8, zorder=0)
            ax[i, j].set_xlabel("t [Myr]")
        ax[i, 0].set_ylabel(f"{name}: {title}\nabundance per H", fontsize=14)
        ax[i, 1].set_ylabel("abundance per H")
        ax[i, 2].set_ylabel("T [K]")
        a = ax[i, 3]
        for sol, col, mk, lab in [("bdf", "0.35", "x", "default BDF"),
                                  ("si", "#D55E00", "o", "semi_implicit")]:
            ws, we = [], []
            for tc in times:
                ref = final(f"{CHECK}/{name}_t{tc}_ref")[1]
                e = errors(final(f"{CHECK}/{name}_t{tc}_{sol}")[1], ref)
                we.append(max(e.pop("eint"), 1e-6))
                ws.append(max(max(e.values()), 1e-6))
            a.loglog([float(v) for v in times], ws, "-", marker=mk, color=col, label=f"{lab}, worst species")
            a.loglog([float(v) for v in times], we, ":", marker=mk, color=col, label=f"{lab}, internal energy")
        a.axhline(0.1, color="0.6", lw=0.8, ls="--")
        a.set_ylim(5e-7, 2)
        a.set_xlabel("check time [Myr]")
        a.set_ylabel("error vs converged BDF")
        a.legend(fontsize=11)
    for j, t in enumerate(["(a) H$_2$, CO, H$^+$", "(b) remaining species",
                           "(c) temperature", "(d) error at check times"]):
        ax[0, j].set_title(t)
    fig.suptitle("GOW17 uniform, 4$^3$: converged Kokkos BDF (rtol = 10$^{-7}$, "
                 "atol = 10$^{-16}$; solid) vs semi_implicit at its defaults (dashed)\n"
                 "default BDF in (d): rtol = 10$^{-3}$, atol = 10$^{-6}$; grey lines: "
                 "check times", fontsize=17)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
