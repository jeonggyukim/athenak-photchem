"""Element streams: the plane-parallel PDR drawn as one band per element.

Each band holds stacked streams, one per species carrying the element, with
thickness equal to that species' share of the element at each depth (GOW17
equilibrium from layouts/slab.py). Species too rare to see at linear scale get
a thin log-abundance strip under the band ("trace threads").
"""
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

from .layouts import slab  # noqa: E402

ABUND = {"C": 1.6e-4, "O": 3.2e-4, "He": 0.1, "Si": 1.7e-6}

# element -> [(species, atoms of the element per molecule, colour)]
STREAMS = {
    "H": [("H+", 1, "#E74C3C"), ("H", 1, "#F5B7B1"), ("H2", 2, "#5DADE2")],
    "He": [("He+", 1, "#E74C3C"), ("He", 1, "#D5D8DC")],
    "C": [("C+", 1, "#E74C3C"), ("C", 1, "#F4D03F"), ("CHx", 1, "#A569BD"),
          ("CO", 1, "#48C9B0"), ("HCO+", 1, "#16A085")],
    "O": [("O+", 1, "#E74C3C"), ("O", 1, "#F8C471"), ("OHx", 1, "#A569BD"),
          ("CO", 1, "#48C9B0"), ("HCO+", 1, "#16A085")],
    "Si": [("Si+", 1, "#E74C3C"), ("Si", 1, "#D5D8DC")],
}
# Trace species shown on a log strip under the band of their element.
TRACES = {
    "H": [("H2+", "#85C1E9"), ("H3+", "#2874A6"), ("e-", "#7F8C8D")],
    "He": [("He+", "#E74C3C")],
    "C": [("CHx", "#A569BD"), ("HCO+", "#16A085")],
    "O": [("OHx", "#A569BD"), ("O+", "#E74C3C")],
    "Si": [],
}
TEX = {"H+": "H$^+$", "H2": "H$_2$", "He+": "He$^+$", "C+": "C$^+$", "CHx": "CH$_x$",
       "HCO+": "HCO$^+$", "O+": "O$^+$", "OHx": "OH$_x$", "Si+": "Si$^+$",
       "H2+": "H$_2^+$", "H3+": "H$_3^+$", "e-": "e$^-$"}
ORDER = ["H", "He", "C", "O", "Si"]
BAND_H, TRACE_H, GAP = 1.0, 0.45, 0.25
LOG_RANGE = (-10.0, -3.0)          # log x shown on trace strips


def partitions(av, x):
    """{element: [(species, fraction array, colour)]}, fractions summing to 1."""
    out = {}
    for elem, parts in STREAMS.items():
        raw = [(s, n * x[s], c) for s, n, c in parts]
        tot = sum(r[1] for r in raw)
        out[elem] = [(s, f / tot, c) for s, f, c in raw]
    return out


def draw(ax=None, av=None, x=None, title=True):
    """Draw the streams on `ax` (new figure if None); return (fig, ax, geometry).

    geometry[elem] = (y_bottom of band, y_top of band, y_bottom of trace strip).
    """
    if av is None:
        av, x = slab.profile()
    fig = None
    if ax is None:
        fig, ax = plt.subplots(figsize=(13.333, 7.5), dpi=200)
    lx = np.log10(av)
    parts = partitions(av, x)
    geom = {}
    y = 0.0
    for elem in reversed(ORDER):
        traces = TRACES[elem]
        y_tr = y
        if traces:
            lo, hi = LOG_RANGE
            ax.add_patch(plt.Rectangle((lx[0], y), lx[-1] - lx[0], TRACE_H,
                                       fc="#FBFCFC", ec="#D5D8DC", lw=0.6))
            for s, c in traces:
                v = np.clip(np.log10(np.maximum(x[s], 1e-30)), lo, hi)
                yy = y + TRACE_H * (v - lo) / (hi - lo)
                ax.plot(lx, yy, color=c, lw=1.6)
                k = int(np.argmax(yy))
                ax.text(lx[k], yy[k] + 0.03, TEX.get(s, s), color=c, fontsize=8,
                        ha="center", va="bottom")
            y += TRACE_H + 0.06
        base = np.full_like(lx, y)
        for s, f, c in parts[elem]:
            top = base + BAND_H * f
            ax.fill_between(lx, base, top, color=c, lw=0, alpha=0.92)
            w = top - base
            k = int(np.argmax(w))
            if w[k] > 0.14:
                kx = min(max(k, 3), len(lx) - 4)   # keep labels off the edges
                ax.text(lx[kx], (base[kx] + top[kx]) / 2, TEX.get(s, s), ha="center",
                        va="center", fontsize=12)
            base = top
        ax.text(lx[0] - 0.06, y + BAND_H / 2, elem, ha="right", va="center",
                fontsize=15, weight="bold")
        if traces:
            ax.text(lx[0] - 0.06, y_tr + TRACE_H / 2, "trace\n(log x)", ha="right",
                    va="center", fontsize=7, color="#7F8C8D")
        geom[elem] = (y, y + BAND_H, y_tr)
        y += BAND_H + GAP

    for name, a in slab.fronts(av, x).items():
        ax.axvline(np.log10(a), color="#2C3E50", ls="--", lw=1)
        ax.text(np.log10(a), y - GAP + 0.08, f"{name}\nA$_V$ = {a:.2f}", ha="center",
                va="bottom", fontsize=9)
    ax.set_xlim(lx[0] - 0.55, lx[-1] + 0.02)
    ax.set_ylim(-0.25, y + 0.55)
    ax.set_xticks(np.log10([0.01, 0.1, 1, 10]))
    ax.set_xticklabels(["0.01", "0.1", "1", "10"])
    ax.set_xlabel("A$_V$ from the illuminated surface (FUV enters from the left)")
    ax.set_yticks([])
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    if title:
        ax.set_title(
            "Plane-parallel PDR, GOW17 at equilibrium ($n_{\\rm H}$ = 10$^3$ cm$^{-3}$, "
            "$\\chi$ = 30): stream thickness = share of each element in each species",
            fontsize=11)
    return fig, ax, geom
