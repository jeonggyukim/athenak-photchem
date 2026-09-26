#!/usr/bin/env python
"""Volume-integrated history of the decaying-turbulence runs, from the dumps, for any
number of runs.

Each run is given with --run TYPE:DIR:BASE[:LABEL], TYPE being ak (AthenaK binary
dumps DIR/bin/BASE.hydro_w.*.bin, timestep from DIR/BASE.hydro.hst) or ap (Athena++
per-MeshBlock VTK DIR/BASE.block*.out1.*.vtk, summed directly so a run can be plotted
while it is still going and before its blocks are joined; timestep from the `cycle=`
lines of DIR/run.out). One run is enough. With several, the first Athena++ run (or the
first run if there is none) is the reference that the final-time differences are
printed against.

Usage:
    plot_turb_history_dumps.py --run ap:turb64_cvode_run:turb64_cvode \
        --run ak:turb64_si_run:turb64_si --run ak:turb64_bdf_run:turb64_bdf:"AthenaK BDF" \
        --out fig.png
The older form (--ak-si-dir/--ak-si-base, --ap-cvode-dir/--ap-cvode-base,
--ak-bdf-dir/--ak-bdf-base, --ak-extra DIR BASE LABEL) still works.
"""
import argparse
import glob
import os
import pathlib
import re
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

GAMMA = 5.0 / 3.0
MYR_PER_CODE = 0.97784
ERG_PER_CODE = 6.882615081124725e31 * 1.0e10  # [mass] (km/s)^2
XC_TOT = 1.6e-4
# Code pressure unit is (1.4 m_H cm^-3)(km/s)^2 and code density is n_H; this turns
# pressure into P/k_B [K cm^-3].
PK_PER_CODE = 1.4 * 1.6735575e-24 * 1.0e10 / 1.380649e-16
IONS = ["He+", "C+", "HCO+", "H+", "H3+", "H2+", "O+", "Si+"]
T_COLD = 300.0  # [K] upper bound of the cold phase
DEFAULT_LABEL = {"ak": "AthenaK semi-implicit", "ap": "Athena++ CVODE"}
STYLES = [dict(ls="--", color="#2E8B57", lw=1.8),
          dict(ls="-", color="#C1440E", lw=4.0, alpha=0.4),
          dict(ls=":", color="#3B5BA9", lw=1.8),
          dict(ls="-.", color="#7A2A06", lw=1.4),
          dict(ls=":", color="#6A3D9A", lw=1.8),
          dict(ls=(0, (5, 1, 1, 1, 1, 1)), color="#B8860B", lw=1.4)]


def integrals(t, dV, rho, v1, v2, v3, press, species):
    m = rho * dV
    M = m.sum()
    # particles per H: H atoms and H2 (1 - x_H2 for atomic + molecular H), He, electrons
    xe = sum(species[s] for s in IONS)
    T = press * PK_PER_CODE / (rho * (1.1 + xe - species["H2"]))
    row = {
        "t": t * MYR_PER_CODE,
        "mass": M,
        "KE": 0.5 * np.sum(m * (v1**2 + v2**2 + v3**2)) * ERG_PER_CODE,
        "Eth": np.sum(press) * dV / (GAMMA - 1.0) * ERG_PER_CODE,
        "fcold": np.sum(m[T < T_COLD]) / M,
        "H2": np.sum(m * 2.0 * species["H2"]) / M,
        "CO": np.sum(m * species["CO"]) / M / XC_TOT,
    }
    for name, x in species.items():
        row["x_" + name] = np.sum(m * x) / M
    return row


def athenak_history(ak_dir, base):
    rows = []
    for f in sorted(glob.glob(f"{ak_dir}/bin/{base}.hydro_w.*.bin")):
        d = bin_convert.read_binary_as_athdf(f)
        dV = np.prod([np.diff(d[k])[0] for k in ("x1f", "x2f", "x3f")])
        g = {k: np.asarray(d[k], dtype=np.float64) for k in
             ("dens", "velx", "vely", "velz", "eint")}
        species = {name: np.asarray(d[f"s_{i:02d}_chem_{name}"], dtype=np.float64)
                   for i, name in enumerate(SPECIES)}
        rows.append(integrals(float(d["Time"]), dV, g["dens"], g["velx"], g["vely"],
                              g["velz"], g["eint"] * (GAMMA - 1.0), species))
    hst = np.loadtxt(glob.glob(f"{ak_dir}/{base}.hydro.hst")[0])
    return rows, hst[:, 0] * MYR_PER_CODE, hst[:, 1]


def vtk_time_and_spacing(fname):
    """Time and uniform cell width from a RECTILINEAR_GRID VTK header."""
    with open(fname, "rb") as fh:
        fh.readline()
        t = float(re.search(r"time=(\S+)", fh.readline().decode("ascii")).group(1))
        line = b""
        while not line.startswith(b"X_COORDINATES"):
            line = fh.readline()
        x = np.frombuffer(fh.read(4 * int(line.split()[1])), dtype=">f4")
    return t, float(x[1] - x[0])


def athenapp_history(ap_dir, base):
    rows = []
    steps = sorted({f.rsplit(".", 2)[-2] for f in
                    glob.glob(f"{ap_dir}/{base}.block*.out1.*.vtk")})
    nblock = len(glob.glob(f"{ap_dir}/{base}.block*.out1.00000.vtk"))
    for n in steps:
        files = sorted(glob.glob(f"{ap_dir}/{base}.block*.out1.{n}.vtk"))
        if len(files) < nblock:
            continue
        t, dx = vtk_time_and_spacing(files[0])
        blocks = [read_athenapp_vtk(f) for f in files]
        keys = ["rho", "vel1", "vel2", "vel3", "press"] + ["r" + s for s in SPECIES]
        cat = {k: np.concatenate([b[k].ravel() for b in blocks]).astype(np.float64)
               for k in keys}
        species = {name: cat["r" + name] for name in SPECIES}
        rows.append(integrals(t, dx**3, cat["rho"], cat["vel1"], cat["vel2"],
                              cat["vel3"], cat["press"], species))
    t_dt = [(float(m.group(1)), float(m.group(2))) for m in
            re.finditer(r"^cycle=\d+ time=(\S+) dt=(\S+)",
                        pathlib.Path(f"{ap_dir}/run.out").read_text(), re.M)]
    t_dt = np.array(t_dt).reshape(-1, 2)
    return rows, t_dt[:, 0] * MYR_PER_CODE, t_dt[:, 1]


def load_run(spec):
    """TYPE:DIR:BASE[:LABEL] -> dict(kind, label, rows, t, dt)."""
    parts = spec.split(":", 3)
    if len(parts) < 3 or parts[0] not in ("ak", "ap"):
        raise SystemExit(f"--run {spec!r}: expected ak|ap:DIR:BASE[:LABEL]")
    kind, directory, base = parts[:3]
    label = parts[3] if len(parts) == 4 else DEFAULT_LABEL[kind]
    rows, t, dt = (athenak_history if kind == "ak" else athenapp_history)(directory, base)
    if not rows:
        raise SystemExit(f"--run {spec!r}: no dumps found")
    print(f"{label}: {len(rows)} dumps to {rows[-1]['t']:.3f} Myr")
    return dict(kind=kind, label=label, rows=rows, t=t, dt=dt)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", action="append", default=[],
                        help="ak|ap:DIR:BASE[:LABEL]; repeat for each run")
    parser.add_argument("--ak-si-dir")
    parser.add_argument("--ak-si-base")
    parser.add_argument("--ak-si-label", default=DEFAULT_LABEL["ak"])
    parser.add_argument("--ap-cvode-dir")
    parser.add_argument("--ap-cvode-base")
    parser.add_argument("--ak-extra", nargs=3, action="append", default=[],
                        metavar=("DIR", "BASE", "LABEL"),
                        help="another AthenaK run; repeatable")
    parser.add_argument("--ak-bdf-dir")
    parser.add_argument("--ak-bdf-base")
    parser.add_argument("--tmax", type=float, help="plot only t <= tmax [Myr]")
    parser.add_argument("--title", default="Decaying turbulence from a shared snapshot")
    parser.add_argument("--out", default="fig_turb_history_dumps.png")
    args = parser.parse_args()

    specs = list(args.run)
    if args.ap_cvode_dir and args.ap_cvode_base:
        specs.append(f"ap:{args.ap_cvode_dir}:{args.ap_cvode_base}")
    if args.ak_si_dir and args.ak_si_base:
        specs.append(f"ak:{args.ak_si_dir}:{args.ak_si_base}:{args.ak_si_label}")
    for d, base, label in args.ak_extra:
        specs.append(f"ak:{d}:{base}:{label}")
    if args.ak_bdf_dir and args.ak_bdf_base:
        specs.append(f"ak:{args.ak_bdf_dir}:{args.ak_bdf_base}:AthenaK Kokkos BDF")
    if not specs:
        parser.error("give at least one --run")

    runs = [load_run(s) for s in specs]
    # Truncating the data, not just the axis, lets the log y-axes rescale to the
    # plotted interval.
    if args.tmax is not None:
        for r in runs:
            keep = r["t"] <= args.tmax
            r["rows"] = [row for row in r["rows"] if row["t"] <= args.tmax]
            r["t"], r["dt"] = r["t"][keep], r["dt"][keep]
    col = lambda rows, k: np.array([r[k] for r in rows])  # noqa: E731
    tmax = max(r["rows"][-1]["t"] for r in runs)

    panels = [
        ("KE", "kinetic energy [erg]", "log"),
        ("Eth", "thermal energy [erg]", "log"),
        ("fcold", rf"cold-gas mass fraction, $T < {T_COLD:.0f}$ K", "linear"),
        ("H2", r"mass-weighted $2x({\rm H_2})$", "log"),
        ("CO", r"mass-weighted $x({\rm CO})/x_{\rm C,tot}$", "log"),
        ("dt", r"$\Delta t$ [Myr]", "log"),
    ]
    fig, axes = plt.subplots(2, 3, figsize=(15, 8), constrained_layout=True)
    for ax, (key, label, scale) in zip(axes.flat, panels):
        for i, r in enumerate(runs):
            style = STYLES[i % len(STYLES)]
            if key == "dt":
                ax.plot(r["t"], r["dt"] * MYR_PER_CODE, **style)
            else:
                ax.plot(col(r["rows"], "t"), col(r["rows"], key), label=r["label"],
                        **style)
        ax.set(xlabel="t [Myr]", ylabel=label, yscale=scale, xlim=(0, tmax))
        ax.grid(alpha=0.25, lw=0.5)
    axes[0, 0].legend(frameon=False, fontsize=9)
    fig.suptitle(f"{args.title}: volume-integrated history")
    fig.savefig(args.out, dpi=140)
    print("wrote", args.out)

    fig, axes = plt.subplots(3, 4, figsize=(18, 11), constrained_layout=True)
    for ax, name in zip(axes.flat, SPECIES):
        for i, r in enumerate(runs):
            ax.plot(col(r["rows"], "t"), col(r["rows"], "x_" + name), label=r["label"],
                    **STYLES[i % len(STYLES)])
        ax.set(xlabel="t [Myr]", ylabel=f"mass-weighted x({name})", yscale="log",
               xlim=(0, tmax))
        ax.grid(alpha=0.25, lw=0.5)
    axes[0, 0].legend(frameon=False, fontsize=8)
    fig.suptitle(f"{args.title}: mass-weighted abundance per H of every integrated "
                 "species")
    out = pathlib.Path(args.out)
    species_out = out.with_name(out.stem + "_species" + out.suffix)
    fig.savefig(species_out, dpi=140)
    print("wrote", species_out)

    if len(runs) > 1:
        ref = next((r for r in runs if r["kind"] == "ap"), runs[0])
        for r in runs:
            if r is ref:
                continue
            for key in ("mass", "KE", "Eth", "H2", "CO"):
                a, p = col(r["rows"], key), col(ref["rows"], key)
                n = min(len(a), len(p))
                rel = (a[n - 1] - p[n - 1]) / p[n - 1]
                print(f"  {r['label']:28s} {key:4s} at t={ref['rows'][n-1]['t']:.3f} Myr:"
                      f" / {ref['label']} - 1 = {rel:+.3e}")


if __name__ == "__main__":
    main()
