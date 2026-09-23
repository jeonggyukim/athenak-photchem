"""Ionization ladder for warm and hot gas.

One column per element, neutral at the bottom and the ion above it, so upward
arrows are ionization and downward arrows recombination. The molecules sit in
a compact block to the right of the ladder.

Columns run Si, He, H, O, C so that He+ -> H+ and H+ <-> O+ join neighbours and
carbon, with the most links to the molecules, sits next to the block. Links
from the ion row to the block pass over the ladder, links from the neutral row
pass under it.
"""
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from . import Layout
from ..reactions import index_of

TITLE = "warm and hot gas: ionization ladder"

_NEU, _ION = 2.5, 7.0
_COL = {"Si": 0.0, "He": 3.0, "H": 6.0, "O": 9.0, "C": 12.0}

_POS = {}
for _el, _x in _COL.items():
    _POS[_el] = (_x, _NEU)
    _POS[_el + "+"] = (_x, _ION)
_POS.update({
    "CHx": (14.2, 4.1), "CO": (15.6, _NEU), "HCO+": (16.6, 6.5),
    "OHx": (15.2, 9.4), "H3+": (18.4, 4.75), "H2+": (20.2, _ION), "H2": (20.2, _NEU),
})

_FRAME = dict(xlim=(-2.6, 29.8), ylim=(-3.2, 12.0), legend_xy=(22.0, 9.6),
              equation_xy=(9.0, -1.8), caption_xy=(9.0, -2.6))


def _background(ax, layout):
    x0 = min(_COL.values()) - 1.6
    ax.add_patch(FancyArrowPatch((x0, _NEU - 0.4), (x0, _ION + 0.4),
                                 arrowstyle="-|>,head_length=8,head_width=4",
                                 lw=2.0, color="#AAB7B8", zorder=0))
    ax.text(x0 - 0.35, (_NEU + _ION) / 2, "ionization", rotation=90,
            ha="center", va="center", fontsize=11, color="#7F8C8D")
    xl, xr = min(_COL.values()) - 1.0, max(x for x, _ in _POS.values()) + 0.9
    yb = 10.2
    ax.add_patch(FancyBboxPatch((xl, yb), xr - xl, 0.9,
                                boxstyle="round,pad=0.02,rounding_size=0.15",
                                fc="#FDEBD0", ec="none", zorder=0))
    ax.text((xl + xr) / 2, yb + 0.45,
            r"above $2\times10^4$ K, NCR takes the He and metal ionization "
            "state from CIE tables",
            ha="center", va="center", fontsize=10.5, color="#6E2C00", zorder=0.5)


_ARCS = {
    # He+ -> H2+ and H2+ -> H+ arc over the ion row.
    index_of("2body", 19): -0.24, index_of("2body", 26): 0.21,
    # O+ <-> OHx both bend upward to clear C+.
    index_of("2body", 29): 0.04, index_of("2body", 25): 0.14,
    # OHx -> O passes right of C+.
    index_of("photo", 3): -0.1, index_of("2body", 21): -0.25,
    # H <-> H2 arc under the neutral row.
    index_of("grain", 0): 0.26, index_of("photo", 4): -0.33,
    index_of("2body", 15): -0.4,
    # NCR's EUV destruction of H2 runs inside the bundle, between UV and collisional.
    index_of("ncr", 1): -0.2,
    # H3+ -> H passes above CHx; CO -> C+ passes right of CHx.
    index_of("2body", 18): 0.25, index_of("2body", 4): 0.5,
    # OHx -> CO passes left of HCO+; H3+ -> OHx passes right of it.
    index_of("2body", 8): 0.08, index_of("2body", 1): 0.35,
}
_LABELS = {index_of("2body", 26): 0.35, index_of("2body", 20): 0.35,
           index_of("2body", 21): 0.86, index_of("2body", 8): 0.3, index_of("2body", 0): 0.25,
           index_of("2body", 5): 0.65, index_of("2body", 29): 0.3}

LADDER = Layout("ladder", dict(_POS), background=_background,
                arc_overrides=_ARCS, label_overrides=_LABELS, **_FRAME)
