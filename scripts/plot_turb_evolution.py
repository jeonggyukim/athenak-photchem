#!/usr/bin/env python
"""Slices, movie and history for the decaying-turbulence runs, for any number of runs.

Each run is given with --run TYPE:DIR:BASE[:LABEL], TYPE being ak (AthenaK binary
dumps in DIR/bin/BASE.hydro_w.*.bin) or ap (Athena++ joined VTK in
DIR/BASE.block0.out1.*.vtk). One run gives single-run slices; several give one row
per run for each group of four quantities. All runs share one colour scale per
quantity, fixed across every frame and every run, so that what changes in the movie
is the gas and not the normalisation. Frames follow the dump times of the first run;
every other run shows its dump nearest in time.

Usage:
    plot_turb_evolution.py --run ap:turb_cvode_run/joined:turb_cvode:"Athena++ CVODE" \
        --run ak:turb_si_run:turb_si:"AthenaK semi-implicit" --outdir figures
The older two-run form (--ak-dir --ak-base --ap-dir --ap-base) still works.
"""
import argparse
import glob
import os
import pathlib
import subprocess
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ATHENAK_ROOT = pathlib.Path(
    os.environ.get("ATHENAK_ROOT", pathlib.Path.home() / "Projects/athenak-chem")
)
sys.path.insert(0, str(ATHENAK_ROOT / "vis/python"))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import bin_convert  # noqa: E402
from compare_ak_ap import read_athenapp_vtk  # noqa: E402

GAMMA = 5.0 / 3.0
MYR_PER_CODE = 0.97784  # pc/(km/s) expressed in Myr
# Code pressure unit is (1.4 m_H cm^-3)(km/s)^2; this converts it to P/k_B.
PK_PER_CODE = 1.4 * 1.6735575e-24 * 1.0e10 / 1.380649e-16
XC_TOT = 1.6e-4
QUANTITIES = [
    ("nH", r"$n_{\rm H}$ [cm$^{-3}$]", "viridis"),
    ("Pk", r"$P/k_{\rm B}$ [cm$^{-3}$ K]", "inferno"),
    ("H2", r"$2x({\rm H_2})$", "cividis"),
    ("CO", r"$x({\rm CO})/x_{\rm C,tot}$", "magma"),
    ("C+", r"$x({\rm C^+})/x_{\rm C,tot}$", "plasma"),
    ("e", r"$x_{\rm e}$", "viridis"),
    ("H+", r"$x({\rm H^+})$", "cividis"),
    ("H2+", r"$x({\rm H_2^+})$", "magma"),
]
IONS = ["He+", "C+", "HCO+", "H+", "H3+", "H2+", "O+", "Si+"]
DEFAULT_LABEL = {"ak": "AthenaK semi-implicit", "ap": "Athena++ CVODE"}
LINESTYLES = ["-", "--", ":", "-."]


def derived(sp, H2, CO):
    """Panel quantities from a {species: abundance} dict."""
    return {
        "H2": 2.0 * H2,
        "CO": CO / XC_TOT,
        "C+": sp["C+"] / XC_TOT,
        "e": sum(sp[s] for s in IONS),
        "H+": sp["H+"],
        "H2+": sp["H2+"],
    }


def athenak_frame(fname):
    d = bin_convert.read_binary_as_athdf(fname)
    nH = np.asarray(d["dens"], dtype=np.float64)
    press = np.asarray(d["eint"], dtype=np.float64) * (GAMMA - 1.0)
    sp = {k.split("_chem_")[1]: np.asarray(v, dtype=np.float64)
          for k, v in d.items() if k.startswith("s_")}
    return {
        "t": float(d["Time"]) * MYR_PER_CODE,
        "nH": nH,
        "Pk": press * PK_PER_CODE,
        **derived(sp, sp["H2"], sp["CO"]),
    }


def athenapp_frame(fname):
    f = read_athenapp_vtk(fname)
    with open(fname, "rb") as fh:
        fh.readline()
        header = fh.readline().decode("ascii", "replace")
    t = np.nan
    for token in header.split():
        if token.startswith("time="):
            t = float(token.split("=")[1]) * MYR_PER_CODE
    return {
        "t": t,
        "nH": f["rho"],
        "Pk": f["press"] * PK_PER_CODE,
        **derived({s: f["r" + s] for s in IONS}, f["rH2"], f["rCO"]),
    }


def load_run(spec):
    """TYPE:DIR:BASE[:LABEL] -> (label, frames)."""
    parts = spec.split(":", 3)
    if len(parts) < 3 or parts[0] not in ("ak", "ap"):
        raise SystemExit(f"--run {spec!r}: expected ak|ap:DIR:BASE[:LABEL]")
    kind, directory, base = parts[:3]
    label = parts[3] if len(parts) == 4 else DEFAULT_LABEL[kind]
    if kind == "ak":
        files = sorted(glob.glob(f"{directory}/bin/{base}.hydro_w.*.bin"))
        frames = [athenak_frame(f) for f in files]
    else:
        files = sorted(glob.glob(f"{directory}/{base}.block0.out1.*.vtk"))
        frames = [athenapp_frame(f) for f in files]
    if not frames:
        raise SystemExit(f"--run {spec!r}: no dumps found")
    print(f"{label}: {len(frames)} frames to {frames[-1]['t']:.3f} Myr")
    return label, frames


def fixed_ranges(all_frames):
    """One colour range per quantity, over every frame of every run."""
    ranges = {}
    for key, _, _ in QUANTITIES:
        # t = 0 holds the uniform r_init abundances, which would set the floor.
        data = np.concatenate([np.asarray(fr[key]).ravel() for fr in all_frames
                               if fr["t"] > 0])
        data = data[np.isfinite(data) & (data > 0)]
        lo, hi = np.percentile(data, [0.5, 99.9])
        ranges[key] = (max(lo, hi * 1e-8), hi)
    return ranges


def history(frames):
    """Mass-weighted means: what a molecular-gas observation would average."""
    out = {k: [] for k in ("t", "H2", "CO", "Pk")}
    for fr in frames:
        w = fr["nH"]
        out["t"].append(fr["t"])
        for key in ("H2", "CO", "Pk"):
            out[key].append(np.sum(w * fr[key]) / np.sum(w))
    return {k: np.array(v) for k, v in out.items()}


def slice_panel(frames, labels, ranges, outfile, suptitle):
    """For each group of four quantities, one row per run."""
    ncol = 4
    nrun = len(frames)
    ngroup = (len(QUANTITIES) + ncol - 1) // ncol
    nrow = ngroup * nrun
    fig, axes = plt.subplots(nrow, ncol, figsize=(17, 4.2 * nrow),
                             constrained_layout=True, squeeze=False)
    for q, (key, label, cmap) in enumerate(QUANTITIES):
        lo, hi = ranges[key]
        norm = mcolors.LogNorm(vmin=lo, vmax=hi)
        col = q % ncol
        for r, fr in enumerate(frames):
            ax = axes[nrun * (q // ncol) + r, col]
            k = fr[key].shape[0] // 2
            im = ax.imshow(np.clip(fr[key][k], lo, hi), origin="lower", cmap=cmap,
                           norm=norm, extent=(-16, 16, -16, 16))
            if r == 0:
                ax.set_title(label)
            fig.colorbar(im, ax=ax, fraction=0.046)
    for row in range(nrow):
        axes[row, 0].set_ylabel(f"{labels[row % nrun]}\ny [pc]")
    for ax in axes[-1]:
        ax.set_xlabel("x [pc]")
    fig.suptitle(suptitle, fontsize=13)
    fig.savefig(outfile, dpi=130)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run", action="append", default=[],
                        help="ak|ap:DIR:BASE[:LABEL]; repeat for each run")
    parser.add_argument("--ak-dir")
    parser.add_argument("--ak-base")
    parser.add_argument("--ap-dir")
    parser.add_argument("--ap-base")
    parser.add_argument("--title", default="Decaying turbulence, GOW17")
    parser.add_argument("--outdir", default="figures")
    args = parser.parse_args()

    specs = list(args.run)
    if args.ap_dir and args.ap_base:
        specs.append(f"ap:{args.ap_dir}:{args.ap_base}")
    if args.ak_dir and args.ak_base:
        specs.append(f"ak:{args.ak_dir}:{args.ak_base}")
    if not specs:
        parser.error("give at least one --run (or --ak-dir/--ak-base)")

    outdir = pathlib.Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    runs = [load_run(s) for s in specs]
    labels = [lab for lab, _ in runs]
    ranges = fixed_ranges([fr for _, frames in runs for fr in frames])
    for key, _, _ in QUANTITIES:
        print(f"  colour range {key}: {ranges[key][0]:.3e} to {ranges[key][1]:.3e}")

    # Frames follow the first run's dumps; the codes choose their own CFL steps,
    # so every other run shows the dump nearest in time.
    ref = runs[0][1]
    times = [np.array([fr["t"] for fr in frames]) for _, frames in runs]
    sets = [[frames[int(np.argmin(np.abs(tt - fr0["t"])))]
             for (_, frames), tt in zip(runs, times)] for fr0 in ref]

    framedir = outdir / "compare_frames"
    framedir.mkdir(parents=True, exist_ok=True)
    for old in framedir.glob("frame_*.png"):
        old.unlink()
    for i, frames in enumerate(sets):
        slice_panel(frames, labels, ranges, framedir / f"frame_{i:04d}.png",
                    f"{args.title} — t = {frames[0]['t']:.3f} Myr")
    movie = outdir / "compare_slices.mp4"
    subprocess.run(
        ["ffmpeg", "-y", "-loglevel", "error", "-framerate", "8",
         "-i", str(framedir / "frame_%04d.png"), "-pix_fmt", "yuv420p",
         "-vf", "scale=1700:-2", str(movie)],
        check=True,
    )
    print("wrote", movie)

    slice_panel(sets[-1], labels, ranges, outdir / "fig_compare_slices_final.png",
                f"{args.title} — t = {sets[-1][0]['t']:.3f} Myr")
    print("wrote", outdir / "fig_compare_slices_final.png")

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    for i, (label, frames) in enumerate(runs):
        h = history(frames)
        for ax, key in zip(axes, ("H2", "CO", "Pk")):
            ax.plot(h["t"], h[key], LINESTYLES[i % len(LINESTYLES)], lw=1.9,
                    color=f"C{i}", label=label)
    for ax, label in zip(axes, [r"mass-weighted $2x({\rm H_2})$",
                                r"mass-weighted $x({\rm CO})/x_{\rm C,tot}$",
                                r"mass-weighted $P/k_{\rm B}$ [cm$^{-3}$ K]"]):
        ax.set(xlabel="t [Myr]", ylabel=label, yscale="log")
        ax.grid(alpha=0.25, lw=0.5)
    axes[0].legend(frameon=False, fontsize=9)
    nx = ref[-1]["nH"].shape[-1]
    fig.suptitle(f"{args.title}, {ref[-1]['t']:.0f} Myr, {nx}$^3$")
    fig.savefig(outdir / "fig_turb_history.png", dpi=140)
    print("wrote", outdir / "fig_turb_history.png")


if __name__ == "__main__":
    main()
