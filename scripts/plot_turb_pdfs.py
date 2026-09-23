#!/usr/bin/env python
"""Mass-weighted PDFs of the decaying-turbulence runs at chosen times.

One figure per time. Each panel overlays Athena++ CVODE and one or two AthenaK
semi-implicit runs, binned identically. Temperature is computed the same way
for every run, from the pressure, n_H, x(H2) and the electron abundance (the
sum of the ions), so differences reflect the state and not two temperature
formulas.

Usage:
    plot_turb_pdfs.py --ak-si-dir turb64_si_run --ak-si-base turb64_si \
        --ap-cvode-dir turb64_cvode_run --ap-cvode-base turb64_cvode \
        [--ak-si2-dir ... --ak-si2-base ... --ak-si2-label ...] \
        --times 0.5 2.0 --out fig_turb_pdfs.png
"""
import argparse
import glob
import os
import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ATHENAK_ROOT = pathlib.Path(
    os.environ.get("ATHENAK_ROOT", pathlib.Path.home() / "Projects/athenak-chem")
)
sys.path.insert(0, str(ATHENAK_ROOT / "vis/python"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bin_convert  # noqa: E402
from compare_ak_ap import SPECIES, read_athenapp_vtk  # noqa: E402

MYR_PER_CODE = 0.97784
PK_PER_CODE = 1.4 * 1.6735575e-24 * 1.0e10 / 1.380649e-16  # P/k in K cm^-3
XHE, XC_TOT = 0.1, 1.6e-4
IONS = ["He+", "C+", "HCO+", "H+", "H3+", "H2+", "O+", "Si+"]

# (key, axis label, log10 range)
QUANTITIES = [
    ("nH", r"$\log_{10} n_{\rm H}$ [cm$^{-3}$]", (-1.5, 3.0)),
    ("Pk", r"$\log_{10} P/k$ [K cm$^{-3}$]", (1.5, 5.5)),
    ("T", r"$\log_{10} T$ [K]", (0.5, 5.0)),
    ("H2", r"$\log_{10} 2x({\rm H_2})$", (-6.0, 0.1)),
    ("CO", r"$\log_{10} x({\rm CO})/x_{\rm C,tot}$", (-9.0, 0.1)),
    ("C+", r"$\log_{10} x({\rm C^+})/x_{\rm C,tot}$", (-6.0, 0.1)),
    ("e", r"$\log_{10} x_e$", (-6.0, -1.0)),
    ("H+", r"$\log_{10} x({\rm H^+})$", (-8.0, -1.0)),
]


def derived(rho, press, x):
    nH = rho  # code density 1 = 1.4 m_H cm^-3 = n_H of 1 cm^-3
    xe = sum(x[s] for s in IONS)
    Pk = press * PK_PER_CODE
    T = Pk / (nH * (1.0 + XHE - x["H2"] + xe))
    return {"w": rho, "nH": nH, "Pk": Pk, "T": T, "H2": 2.0 * x["H2"],
            "CO": x["CO"] / XC_TOT, "C+": x["C+"] / XC_TOT, "e": xe, "H+": x["H+"]}


def athenak_state(run_dir, base, t_myr):
    files = sorted(glob.glob(f"{run_dir}/bin/{base}.hydro_w.*.bin"))
    times = [float(bin_convert.read_binary_as_athdf(f)["Time"]) * MYR_PER_CODE
             for f in files]
    f = files[int(np.argmin(np.abs(np.array(times) - t_myr)))]
    d = bin_convert.read_binary_as_athdf(f)
    x = {s: np.asarray(d[f"s_{i:02d}_chem_{s}"], dtype=float).ravel()
         for i, s in enumerate(SPECIES)}
    rho = np.asarray(d["dens"], dtype=float).ravel()
    press = np.asarray(d["eint"], dtype=float).ravel() * (5.0 / 3.0 - 1.0)
    return float(d["Time"]) * MYR_PER_CODE, derived(rho, press, x)


def athenapp_state(run_dir, base, t_myr):
    steps = sorted({f.rsplit(".", 2)[-2] for f in
                    glob.glob(f"{run_dir}/{base}.block0.out1.*.vtk")})
    times = []
    for n in steps:
        with open(f"{run_dir}/{base}.block0.out1.{n}.vtk", "rb") as fh:
            fh.readline()
            line = fh.readline().decode("ascii")
        times.append(float(line.split("time=")[1].split()[0]) * MYR_PER_CODE)
    n = steps[int(np.argmin(np.abs(np.array(times) - t_myr)))]
    blocks = [read_athenapp_vtk(f) for f in
              sorted(glob.glob(f"{run_dir}/{base}.block*.out1.{n}.vtk"))]
    cat = lambda k: np.concatenate([b[k].ravel() for b in blocks]).astype(float)  # noqa: E731
    x = {s: cat("r" + s) for s in SPECIES}
    return times[steps.index(n)], derived(cat("rho"), cat("press"), x)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ak-si-dir", required=True)
    parser.add_argument("--ak-si-base", required=True)
    parser.add_argument("--ap-cvode-dir", required=True)
    parser.add_argument("--ap-cvode-base", required=True)
    parser.add_argument("--ak-si2-dir")
    parser.add_argument("--ak-si2-base")
    parser.add_argument("--ak-si2-label", default="AthenaK semi-implicit, variant")
    parser.add_argument("--times", type=float, nargs="+", default=[0.5, 2.0],
                        help="times in Myr; the nearest dump is used")
    parser.add_argument("--fixed-range", action="store_true",
                        help="use the hand-set ranges, comparable across figures, "
                             "instead of the data's min and max")
    parser.add_argument("--out", default="fig_turb_pdfs.png")
    args = parser.parse_args()

    out = pathlib.Path(args.out)
    for t in args.times:
        runs = [("Athena++ CVODE", athenapp_state(args.ap_cvode_dir,
                                                  args.ap_cvode_base, t),
                 dict(ls="--", color="#2E8B57", lw=1.8)),
                ("AthenaK semi-implicit", athenak_state(args.ak_si_dir,
                                                        args.ak_si_base, t),
                 dict(ls="-", color="#C1440E", lw=4.0, alpha=0.4))]
        if args.ak_si2_dir:
            runs.append((args.ak_si2_label,
                         athenak_state(args.ak_si2_dir, args.ak_si2_base, t),
                         dict(ls="-.", color="#7A2A06", lw=1.4)))

        fig, axes = plt.subplots(2, 4, figsize=(18, 8), constrained_layout=True)
        for ax, (key, label, (lo, hi)) in zip(axes.flat, QUANTITIES):
            logs = [np.log10(np.clip(s[key], 1e-300, None)) for _, (_, s), _ in runs]
            if not args.fixed_range:
                lo = min(v.min() for v in logs)
                hi = max(v.max() for v in logs)
                pad = max(0.02 * (hi - lo), 0.5 * (0.2 - (hi - lo)), 0.0)
                lo, hi = lo - pad, hi + pad
            bins = np.linspace(lo, hi, 81)
            for (name, (t_run, s), style), v in zip(runs, logs):
                h, e = np.histogram(v, bins=bins, weights=s["w"], density=True)
                ax.stairs(h, e, label=f"{name} (t={t_run:.2f} Myr)", **style)
            ax.set(xlabel=label, ylabel="mass-weighted PDF", yscale="log")
            ax.grid(alpha=0.25, lw=0.5)
        axes[0, 0].legend(frameon=False, fontsize=8)
        fig.suptitle(f"Decaying turbulence from a shared snapshot: "
                     f"mass-weighted PDFs at t = {t:.2f} Myr")
        fname = out.with_name(f"{out.stem}_t{t:.2f}{out.suffix}")
        fig.savefig(fname, dpi=140)
        plt.close(fig)
        print("wrote", fname)


if __name__ == "__main__":
    main()
