"""Movie: carbon atoms settling into the PDR layers.

Every depth column of the slab holds the same number of carbon atoms, all C+
at t = 0. Each atom hops between C+, C, CHx, CO and HCO+ with the rates GOW17
gives at its depth (gow17net.swarm). The layers emerge on their own: the
surface settles within years, the C/CO front only after ~1e5 yr.

    python movie/swarm_movie.py            # writes figures/gow17net/movie/swarm.mp4
"""
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

from gow17net import swarm  # noqa: E402
from gow17net.layouts import slab  # noqa: E402

COLORS = {"C+": "#E74C3C", "C": "#F4D03F", "CHx": "#A569BD", "CO": "#48C9B0",
          "HCO+": "#16A085"}
TEX = {"C+": "C$^+$", "C": "C", "CHx": "CH$_x$", "CO": "CO", "HCO+": "HCO$^+$"}
BG, FG = "#101418", "#ECF0F1"


def main():
    n_atoms = 70
    times = np.logspace(0, 6.3, 240)
    av, t, states = swarm.simulate(n_atoms=n_atoms, times_yr=times)
    lx = np.log10(av)
    rng = np.random.default_rng(3)
    jx = rng.uniform(-0.028, 0.028, (len(av), n_atoms))
    jy = rng.uniform(0.02, 0.98, (len(av), n_atoms))
    _, xeq = slab.profile()
    xc = sum(xeq[s] for s in swarm.CARBON)
    fronts = slab.fronts(*slab.profile())
    cols = np.array([matplotlib.colors.to_rgb(COLORS[s]) for s in swarm.CARBON])

    out = ROOT / "figures" / "gow17net" / "movie" / "swarm.mp4"
    out.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as tmp:
        for f in range(len(t)):
            fig = plt.figure(figsize=(16, 9), dpi=120, facecolor=BG)
            ax = fig.add_axes([0.07, 0.36, 0.86, 0.46], facecolor=BG)
            bx = fig.add_axes([0.07, 0.10, 0.86, 0.2], facecolor=BG)
            st = states[f]
            ax.scatter((lx[:, None] + jx).ravel(), jy.ravel(), s=9,
                       c=cols[st.ravel()], linewidths=0)
            for name, a in fronts.items():
                for a_ in (ax, bx):
                    a_.axvline(np.log10(a), color="#7F8C8D", ls="--", lw=1)
                ax.text(np.log10(a), 1.02, f"{name}\nA$_V$ = {a:.2f}", color="#AEB6BF",
                        ha="center", va="bottom", fontsize=11)
            # Stacked fractions now (filled) against equilibrium (lines).
            frac = np.stack([(st == i).mean(axis=1) for i in range(len(swarm.CARBON))])
            base = np.zeros(len(av))
            for i, s in enumerate(swarm.CARBON):
                bx.fill_between(lx, base, base + frac[i], color=COLORS[s], lw=0)
                base = base + frac[i]
            base = np.zeros(len(av))
            for s in swarm.CARBON:
                base = base + xeq[s] / xc
                bx.plot(lx, base, color=FG, lw=0.8, ls=":")
            for a_ in (ax, bx):
                a_.set_xlim(lx[0] - 0.05, lx[-1] + 0.05)
                a_.tick_params(colors=FG)
                for sp in a_.spines.values():
                    sp.set_visible(False)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_xticks([])
            bx.set_ylim(0, 1)
            bx.set_xticks(np.log10([0.01, 0.1, 1, 10]))
            bx.set_xticklabels(["0.01", "0.1", "1", "10"], color=FG)
            bx.set_xlabel("A$_V$ from the illuminated surface", color=FG, fontsize=13)
            bx.set_ylabel("carbon\nfraction", color=FG, fontsize=11)
            bx.text(1.0, 1.04, "dotted: equilibrium", transform=bx.transAxes, ha="right",
                    color="#AEB6BF", fontsize=10)
            fig.text(0.07, 0.95, "Carbon atoms settling into a plane-parallel PDR",
                     color=FG, fontsize=20, va="center")
            fig.text(0.93, 0.95, f"t = {t[f]:.3g} yr", color=FG, fontsize=20,
                     ha="right", va="center", family="monospace")
            fig.text(0.07, 0.91, "each dot is a carbon atom; all start as C$^+$ and hop "
                     "between species at the GOW17 rates of their depth "
                     "(n$_{\\rm H}$ = 10$^3$ cm$^{-3}$, $\\chi$ = 30)", color="#AEB6BF",
                     fontsize=12, va="center")
            handles = [plt.Line2D([], [], marker="o", ls="", color=COLORS[s], ms=9,
                                  label=TEX[s]) for s in swarm.CARBON]
            fig.legend(handles=handles, loc="lower center", ncol=5, frameon=False,
                       labelcolor=FG, fontsize=13, bbox_to_anchor=(0.5, 0.0))
            fig.text(0.02, 0.62, "FUV $\\rightarrow$", color="#AF7AC5", fontsize=16,
                     va="center")
            fig.savefig(f"{tmp}/f{f:04d}.png", facecolor=BG)
            plt.close(fig)
        subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", "24",
                        "-i", f"{tmp}/f%04d.png", "-vf", "scale=1920:1080",
                        "-c:v", "libx264", "-pix_fmt", "yuv420p", str(out)], check=True)
    print("wrote", out)


if __name__ == "__main__":
    main()
