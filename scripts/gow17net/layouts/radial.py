"""The classical spherical PDR: a uniform cloud lit by an isotropic external
FUV field. Radius stands for depth: the C+ zone at the surface, the C zone
inside it, the CO zone in the core (fronts from pdr_profile.py at
n_H = 1e3 cm^-3, chi = 1: C+/C at A_V = 0.56, C/CO at A_V = 1.11).

Angle groups the element families: hydrogen left, helium upper left, silicon
lower left, carbon upper right, oxygen lower right, CO and HCO+ facing right,
between carbon and oxygen.

Shell assignment follows the zone each species characterises. O, He and Si are
present at every depth, and O+ peaks deep in the equilibrium profile (cosmic
rays make it there), so their shells are a choice, not a computed peak.
"""
import math

from matplotlib.patches import Circle, FancyArrowPatch

from . import Layout
from ..reactions import index_of

R_CORE, R_MID, R_OUT = 2.4, 5.2, 8.0
ZONE_EDGES = (3.6, 6.6, 9.3)     # CO | C | C+ zone boundaries
LABEL_DEG = -98                  # zone labels sit in the gap between O and Si
FRONTS_AV = {"C/CO": 1.11, "C+/C": 0.56}

_POLAR = {
    # core: CO zone
    "H2": (R_CORE, 150), "H3+": (R_CORE, 210), "CO": (R_CORE, 30), "HCO+": (R_CORE, -30),
    # C zone
    "H2+": (R_MID, 180), "CHx": (R_MID, 45), "C": (R_MID, 78), "OHx": (R_MID, -45),
    "O": (R_MID, -72), "He": (R_MID, 128), "Si": (R_MID, 212),
    # C+ zone, the surface
    "H": (R_OUT, 162), "H+": (R_OUT, 198), "C+": (R_OUT, 58), "O+": (R_OUT, -58),
    "He+": (R_OUT, 120), "Si+": (R_OUT, 228),
}


def _xy(r, deg):
    a = math.radians(deg)
    return (round(r * math.cos(a), 3), round(r * math.sin(a), 3))


def _background(ax, layout):
    shades = ("#F4F6F7", "#E5E8E8", "#D5D8DC")   # C+ zone, C zone, CO zone
    for r, c in zip(reversed(ZONE_EDGES), shades):
        ax.add_patch(Circle((0, 0), r, fc=c, ec="none", zorder=0))
    for r in ZONE_EDGES[:2]:
        ax.add_patch(Circle((0, 0), r, fill=False, ec="#7F8C8D", ls=(0, (4, 3)),
                            lw=1.0, zorder=0.5))
    a = math.radians(LABEL_DEG)
    inner = (0.0,) + ZONE_EDGES[:-1]
    labels = ("CO zone", "C zone\n" + f"C/CO front A$_V$ = {FRONTS_AV['C/CO']}",
              "C$^+$ zone\n" + f"C$^+$/C front A$_V$ = {FRONTS_AV['C+/C']}")
    for r0, r1, text in zip(inner, ZONE_EDGES, labels):
        r = (r0 + r1) / 2 if r0 else 1.1
        ax.text(r * math.cos(a), r * math.sin(a), text, ha="center", va="center",
                fontsize=8, color="#7B7D7D", style="italic", zorder=0.6)
    # FUV enters from outside; cosmic rays reach the core.
    for deg in range(0, 360, 45):
        a = math.radians(deg + 22.5)
        p0 = (10.5 * math.cos(a), 10.5 * math.sin(a))
        p1 = (9.45 * math.cos(a), 9.45 * math.sin(a))
        ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>,head_length=5,head_width=3",
                                     color="#7D3C98", lw=1.4, zorder=0.5))
    ax.text(-11.2, 10.9, "FUV from outside\ncosmic rays reach the core", fontsize=9,
            color="#7D3C98", va="top")


RADIAL = Layout(
    "radial",
    {name: _xy(r, deg) for name, (r, deg) in _POLAR.items()},
    xlim=(-11.4, 21.0), ylim=(-12.9, 11.1), box_w=1.9, box_h=1.0,
    legend_xy=(12.0, 10.2), equation_xy=(4.8, -11.3), caption_xy=(4.8, -12.3),
    background=_background,
    arc_overrides={index_of("2body", 28): -0.05, index_of("2body", 27): 0.3},
)
