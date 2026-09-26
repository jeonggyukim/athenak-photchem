#!/usr/bin/env python
"""History of the decaying-turbulence runs from the codes' .hst files, for any number
of runs.

AthenaK and Athena++ write the same history header format (`[n]=name`), with volume
integrals over the box: mass, kinetic energy per direction and total energy. Thermal
energy is total minus kinetic, since neither run has magnetic fields.

Usage:
    plot_turb_hst.py --hst ap.hst:"Athena++ CVODE" --hst ak.hydro.hst:"AthenaK" \
        --out fig_turb_hst.png
Each --hst is FILE[:LABEL]; one is enough. The older --ak FILE --ap FILE still works.
"""
import argparse
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

MYR_PER_CODE = 0.97784  # pc/(km/s) expressed in Myr


def read_hst(fname):
    with open(fname) as f:
        for line in f:
            if "[1]=" in line:
                names = [m.group(1) for m in re.finditer(r"\[\d+\]=(\S+)", line)]
                break
    data = np.loadtxt(fname, ndmin=2)
    cols = {name: data[:, i] for i, name in enumerate(names)}
    ke = cols["1-KE"] + cols["2-KE"] + cols["3-KE"]
    return {
        "t": cols["time"] * MYR_PER_CODE,
        "dt": cols["dt"] * MYR_PER_CODE,
        "mass": cols["mass"],
        "KE": ke,
        "TE": cols["tot-E"] - ke,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--hst", action="append", default=[],
                        help="FILE[:LABEL]; repeat for each run")
    parser.add_argument("--ak", help="AthenaK .hydro.hst")
    parser.add_argument("--ap", help="Athena++ .hst")
    parser.add_argument("--out", default="fig_turb_hst.png")
    args = parser.parse_args()

    specs = list(args.hst)
    if args.ap:
        specs.append(f"{args.ap}:Athena++ CVODE")
    if args.ak:
        specs.append(f"{args.ak}:AthenaK semi-implicit")
    if not specs:
        parser.error("give at least one --hst")
    runs = []
    for i, spec in enumerate(specs):
        fname, _, label = spec.partition(":")
        runs.append((label or fname, read_hst(fname), f"C{i}"))

    panels = [("KE", "kinetic energy [code]"),
              ("TE", "thermal energy [code]"),
              ("dt", "dt [Myr]"),
              ("mass", "mass [code]")]
    fig, axes = plt.subplots(2, 2, figsize=(10, 7), sharex=True)
    for ax, (key, label) in zip(axes.flat, panels):
        for name, h, c in runs:
            ax.plot(h["t"], h[key], color=c, label=name)
        ax.set_ylabel(label)
        if key in ("KE", "TE"):
            ax.set_yscale("log")
    for ax in axes[1]:
        ax.set_xlabel("t [Myr]")
    axes[0, 0].legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=140)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
