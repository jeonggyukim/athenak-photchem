"""Two figures on why the slab PDR is layered, from the rates of pdr_rates.py.

fig_rate_crossings.png: at the H/H2, C+/C and C/CO fronts, the formation and
destruction rates per H nucleus of H2, C+ and CO, reaction by reaction.
fig_dominant_process.png: for each key species, the process class with the
largest share of its destruction rate at each A_V.

Run from scripts/:  python -m gow17net.pdr_why
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from . import pdr_rates as P  # noqa: E402
from .layouts import slab  # noqa: E402
from .pdr_profile import OUTDIR  # noqa: E402
from .style import PROCESS_STYLE  # noqa: E402

FIGDIR = OUTDIR / "check" / "pdr_why"

COLOR = {"FUV": PROCESS_STYLE["uv"]["color"],
         "CR": PROCESS_STYLE["cr"]["color"],
         "grain": PROCESS_STYLE["grain"]["color"],
         "recombination": "#17A589",
         "chemistry": PROCESS_STYLE["chem"]["color"],
         "collisional": PROCESS_STYLE["coll"]["color"]}
CLASS_LABEL = {"FUV": "FUV photoreaction", "CR": "cosmic rays (incl. CR-induced photons)",
               "grain": "grain surface", "recombination": "radiative / dissociative recombination",
               "chemistry": "ion-molecule / neutral-neutral", "collisional": "collisional"}
GREY = "#566573"
LS = {"form": ["-", (0, (1, 1))],
      "dest": [(0, (5, 2)), (0, (5, 1.5, 1, 1.5)), (0, (2, 1))]}

# (label, DESPOTIC keys, role, label A_V, label log10 y)
PANELS = [
    ("H/H$_2$", "H2", "H$_2$", [
        ("H + H (grain) → H$_2$", [("tb", 31)], "form", 0.012, -13.2),
        ("H$_2$ + FUV → 2H", [("ph", 5)], "dest", 0.012, -14.05),
        ("H$_2$ + CR → H$_2^+$", [("cr", "H2")], "dest", 1.6, -15.2),
        ("H$_2^+$ + H$_2$ → H$_3^+$", [("tb", 18)], "dest", 2.2, -16.05),
        ("C$^+$ + H$_2$", [("tb", 8), ("tb", 9)], "dest", 0.2, -16.6),
        ("He$^+$ + H$_2$", [("tb", 6), ("tb", 24)], "dest", 2.2, -17.2),
    ]),
    ("C$^+$/C", "C+", "C$^+$", [
        ("C + FUV → C$^+$", [("ph", 0)], "form", 0.012, -14.4),
        ("C + CR → C$^+$", [("cr", 1), ("cr", 3)], "form", 0.03, -19.25),
        ("He$^+$ + CO", [("tb", 7)], "form", 2.0, -16.4),
        ("C$^+$ + e (radiative)", [("tb", 16)], "dest", 0.012, -15.85),
        ("C$^+$ + e (grain)", [("tb", 33)], "dest", 0.012, -15.45),
        ("C$^+$ + H$_2$", [("tb", 8), ("tb", 9)], "dest", 0.2, -15.7),
        ("C$^+$ + OH$_x$ → HCO$^+$", [("tb", 10)], "dest", 0.03, -18.1),
    ]),
    ("C/CO", "CO", "CO", [
        ("CH$_x$ + O → CO", [("tb", 11)], "form", 0.1, -20.65),
        ("OH$_x$ + C → CO", [("tb", 12)], "form", 2.9, -17.15),
        ("HCO$^+$ + e → CO", [("tb", 17)], "form", 0.012, -17.3),
        ("CO + FUV → C + O", [("ph", 2)], "dest", 0.012, -17.75),
        ("CO + CR photons →\nC + O", [("cr", 4)], "dest", 2.9, -18.2),
        ("He$^+$ + CO → C$^+$", [("tb", 7)], "dest", 2.9, -16.45),
        ("H$_3^+$ + CO → HCO$^+$", [("tb", 5)], "dest", 2.9, -15.7),
    ]),
]


def _fronts():
    av, x = slab.profile()
    return slab.fronts(av, x)


def _log_crossing(av, a, b):
    r = np.log(np.maximum(a, 1e-300) / np.maximum(b, 1e-300))
    idx = np.nonzero(np.sign(r[:-1]) != np.sign(r[1:]))[0]
    out = []
    for i in idx:
        la, lb = np.log10(av[i]), np.log10(av[i + 1])
        out.append(10 ** (la + (lb - la) * r[i] / (r[i] - r[i + 1])))
    return out


def fig_rate_crossings(path=FIGDIR / "fig_rate_crossings.png"):
    av, rates = P.rates_vs_av()
    fr = _fronts()
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.8), sharex=True)
    for ax, (front, sp, sp_tex, curves) in zip(axes, PANELS):
        total = sum(P.by_reaction(sp, -1).values())
        ax.plot(av, total, color="#D5D8DC", lw=7, solid_capstyle="butt", zorder=1)
        seen = {}
        lines = {}
        for label, keys, role, av_lab, ly in curves:
            r = sum(rates[k] for k in keys)
            cls = P.process_class(keys[0])
            n = seen.get((cls, role), 0)
            seen[(cls, role)] = n + 1
            ls = LS[role][n]
            ax.plot(av, r, color=COLOR[cls], ls=ls, lw=2.0, zorder=3)
            lines[label] = r
            ax.text(av_lab, 10 ** ly, label, color=COLOR[cls], fontsize=9.5, ha="left",
                    va="center", zorder=5,
                    bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.5))
        ymax = total.max() * 8
        ax.set_ylim(ymax * 1e-6, ymax)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(av[0], av[-1])
        a_f = fr[front]
        ax.axvline(a_f, color=GREY, ls="--", lw=1.2, zorder=2)
        right = front == "C/CO"
        ax.text(a_f * (0.94 if right else 1.06), ymax * 0.6,
                f"{front} front\n$A_V$ = {a_f:.2f}", color=GREY, fontsize=10,
                va="top", ha="right" if right else "left")
        ax.set_title(f"{sp_tex}: formation (solid) vs destruction (dashed)",
                     fontsize=12)
        ax.set_xlabel("$A_V$ [mag]", fontsize=11)
        ax.tick_params(labelsize=10)
        ax.grid(alpha=0.25, which="major")
        _mark_switch(ax, av, curves, lines, sp)
    axes[0].set_ylabel("rate per H nucleus [s$^{-1}$]", fontsize=11)
    handles = [Patch(color=COLOR[c], label=CLASS_LABEL[c])
               for c in ("FUV", "CR", "grain", "recombination", "chemistry")]
    handles += [Line2D([], [], color="k", lw=2, label="formation"),
                Line2D([], [], color="k", lw=2, ls=(0, (5, 2)), label="destruction"),
                Line2D([], [], color="#D5D8DC", lw=7,
                       label="total (formation = destruction at equilibrium)"),
                Line2D([], [], marker="o", ls="none", mfc="white", mec="k",
                       label="leading destruction channel changes")]
    fig.legend(handles=handles, loc="lower center", ncol=5, fontsize=10,
               frameon=False, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Why the fronts sit where they do: rates of the reactions that "
                 "make and destroy H$_2$, C$^+$ and CO   "
                 rf"($n_{{\rm H}} = 10^3$ cm$^{{-3}}$, $\chi = {slab.CHI:g}$, "
                 r"$\xi_{\rm CR} = 2\times10^{-16}$ s$^{-1}$, "
                 rf"$T = {slab.T:.0f}$ K)", fontsize=13)
    fig.tight_layout(rect=(0, 0.09, 1, 0.97))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


def _mark_switch(ax, av, curves, lines, sp):
    """Circle each A_V where the largest destruction curve changes."""
    dest = [(lab, lines[lab]) for lab, _, role, _, _ in curves if role == "dest"]
    arr = np.array([r for _, r in dest])
    lead = np.argmax(arr, axis=0)
    for j in np.nonzero(lead[1:] != lead[:-1])[0]:
        a, b = arr[lead[j]], arr[lead[j + 1]]
        for c in _log_crossing(av[j:j + 2], a[j:j + 2], b[j:j + 2]):
            y = 10 ** np.interp(np.log10(c), np.log10(av), np.log10(a))
            ax.plot(c, y, "o", mfc="white", mec="k", ms=8, mew=1.4, zorder=6)


STRIPS = [("H2", "H$_2$"), ("H+", "H$^+$"), ("C+", "C$^+$"), ("C", "C"),
          ("CO", "CO"), ("e-", "e$^-$")]
SHORT = {"FUV": "FUV", "CR": "CR", "grain": "grain", "recombination": "recomb.",
         "chemistry": "chem.", "collisional": "coll."}


def _pretty(key):
    s = P.label(key).replace("->", "→")
    for a, b in (("H2+", "H$_2^+$"), ("H3+", "H$_3^+$"), ("HCO+", "HCO$^+$"),
                 ("He+", "He$^+$"), ("C+", "C$^+$"), ("H+", "H$^+$"), ("O+", "O$^+$"),
                 ("H2", "H$_2$"), ("CHx", "CH$_x$"), ("OHx", "OH$_x$"), ("e-", "e")):
        s = s.replace(a, b)
    return s


def fig_dominant_process(path=FIGDIR / "fig_dominant_process.png"):
    av, _ = P.rates_vs_av()
    fr = _fronts()
    la = np.log10(av)
    fine = np.linspace(la[0], la[-1], 1200)
    fig, ax = plt.subplots(figsize=(15, 7.4))
    h, gap = 0.78, 0.22
    used = set()
    for row, (sp, tex) in enumerate(STRIPS):
        y0 = (len(STRIPS) - 1 - row) * (h + gap)
        dest = P.destruction_by_process(sp)
        byr = P.by_reaction(sp, -1)
        classes = list(dest)
        lg = np.array([np.interp(fine, la, np.log10(np.maximum(dest[c], 1e-300)))
                       for c in classes])
        tot = np.log10(np.sum(10 ** lg, axis=0))
        share = 10 ** (lg - tot)
        lead = np.argmax(share, axis=0)
        edges = np.concatenate([[0], np.nonzero(lead[1:] != lead[:-1])[0] + 1,
                                [len(fine)]])
        for i0, i1 in zip(edges[:-1], edges[1:]):
            c = classes[lead[i0]]
            x0, x1 = 10 ** fine[i0], 10 ** fine[min(i1, len(fine) - 1)]
            ax.fill_between([x0, x1], y0, y0 + h, color=COLOR[c], alpha=1.0,
                            lw=0, zorder=1)
            width = fine[min(i1, len(fine) - 1)] - fine[i0]
            im = (i0 + i1 - 1) // 2
            s = share[lead[i0], im]
            amid = 10 ** fine[im]
            seg = share[lead[i0], i0:i1]
            lo, hi = 100 * seg.min(), 100 * seg.max()
            pct = f"{lo:.0f}%" if hi - lo < 5 else f"{lo:.0f}-{hi:.0f}%"
            txt = f"{SHORT[c]} {pct}"
            others = [k for k in range(len(classes)) if k != lead[i0]]
            second = max(others, key=lambda k: share[k, im])
            if share[second, im] > 0.3:
                txt += f" (then {SHORT[classes[second]]} {100 * share[second, im]:.0f}%)"
            used.add(c)
            if width > 0.55:
                keys = [k for k in byr if P.process_class(k) == c]
                kbest = max(keys, key=lambda k: np.interp(fine[im], la, byr[k]))
                txt += "\n" + _pretty(kbest)
            fs = 10 if width > 0.3 else 8.5
            if width > 0.12:
                ax.text(amid, y0 + h / 2, txt, ha="center", va="center",
                        fontsize=fs, color="white", zorder=5, linespacing=1.15,
                        bbox=dict(fc=COLOR[c], ec="none", pad=1.5))
        ax.text(10 ** (la[0] - 0.04), y0 + h / 2, tex, ha="right", va="center",
                fontsize=14)
    ytop = len(STRIPS) * (h + gap) - gap
    for name, a in fr.items():
        ax.plot([a, a], [-0.05, ytop + 0.05], ls="--", lw=1.4, color="#1B2631",
                zorder=4)
        ax.text(a, ytop + 0.08, f"{name}\n$A_V$ = {a:.2f}", ha="center",
                va="bottom", fontsize=10, color="#1B2631")
    ax.set_xscale("log")
    ax.set_xlim(av[0], av[-1])
    ax.set_ylim(-0.1, ytop + 0.55)
    ax.set_yticks([])
    for s in ("left", "right", "top"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("$A_V$ [mag], depth from the illuminated surface", fontsize=12)
    ax.tick_params(labelsize=11)
    present = [c for c in COLOR if c in used]
    ax.legend(handles=[Patch(color=COLOR[c], alpha=1.0, label=CLASS_LABEL[c])
                       for c in present], loc="upper center",
              bbox_to_anchor=(0.5, -0.13), ncol=3, fontsize=10.5, frameon=False)
    ax.set_title("Which process controls each species at each depth\n"
                 "colour: process class with the largest share of the species' "
                 "destruction rate; text: that share and its leading reaction",
                 fontsize=13, pad=12)
    fig.text(0.5, 0.012,
             "At depth the chemistry share of H$_2$ destruction is H$_2^+$ + H$_2$, "
             "which follows every cosmic-ray ionization of H$_2$.\nCO lost to "
             "H$_3^+$ + CO returns through HCO$^+$ + e; He$^+$ + CO, also "
             "CR-driven, is the net carbon sink of CO.",
             ha="center", va="bottom", fontsize=10, color=GREY)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


if __name__ == "__main__":
    print(fig_rate_crossings())
    print(fig_dominant_process())
