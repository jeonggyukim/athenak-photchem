#!/usr/bin/env python
"""Diagram of the GOW17 chemical network as implemented in AthenaK.

Reactions follow GOW17Network::CDRates in src/chemistry/network/gow17.hpp,
numbered as there within each group (cr, 2-body, photo, grain). Boxes are
species: integrated ions and neutrals are filled, ghost species (derived from
conservation laws, not integrated) are dashed. The number on each integrated
species is its position in the block Gauss-Seidel update; the two species with
thick borders limit the substep; the shaded pair is solved as one coupled block.
Electrons are a ghost species too and appear as the "e" label on recombination
arrows.

Usage:
    plot_gow17_network.py --out fig_gow17_network.png
"""
import argparse
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

POS = {
    "He": (7.4, 8.0), "He+": (7.4, 5.5),
    "H": (1.0, 8.0), "H2": (4.2, 8.0), "H+": (1.0, 3.0),
    "H2+": (4.2, 5.5), "H3+": (4.2, 3.0),
    "C": (10.6, 8.0), "C+": (13.8, 8.0), "CO": (17.0, 8.0),
    "CHx": (11.0, 5.5), "HCO+": (17.0, 5.0),
    "O": (13.8, 0.6), "OHx": (13.8, 3.0), "O+": (10.6, 0.6),
    "Si": (-1.8, 8.0), "Si+": (-1.8, 5.5),
}
TEX = {"He+": "He$^+$", "H2": "H$_2$", "H+": "H$^+$", "H2+": "H$_2^+$",
       "H3+": "H$_3^+$", "C+": "C$^+$", "CHx": "CH$_x$", "HCO+": "HCO$^+$",
       "OHx": "OH$_x$", "O+": "O$^+$", "Si+": "Si$^+$"}
GHOSTS = {"H", "He", "C", "O", "Si"}
IONS = {"He+", "H+", "H2+", "H3+", "C+", "HCO+", "O+", "Si+"}
GS_ORDER = ["He+", "Si+", "H2+", "H3+", "H+", "O+", "C+", "CHx", "OHx",
            "CO", "HCO+", "H2"]
LIMITERS = {"H2", "CO"}
# TIGRESS-NCR (Kim+23) carries a subset, each species closed differently: H2
# and H+ as ODEs, H, C and O as ghosts from the element budgets, C+ and CO from
# equilibrium formulas, O+ as a fixed ratio to H+. Electrons are a charge sum.
NCR = {"H2": "ODE", "H+": "ODE", "H": "ghost", "C": "ghost", "O": "ghost",
       "C+": "equilibrium", "CO": "equilibrium", "O+": "equilibrium"}
NCR_COLOR = {"ODE": "#F5B041", "ghost": "#AAB7B8", "equilibrium": "#48C9B0"}
NCR_LEGEND = {"ODE": "NCR: non-equilibrium", "ghost": "NCR: ghost (element budget)",
              "equilibrium": "NCR: equilibrium or fixed ratio"}

STYLE = {
    "cr": dict(color="#C0392B", ls="-", lw=1.7),
    "crphoto": dict(color="#C0392B", ls=(0, (4, 2)), lw=1.7),
    "uv": dict(color="#7D3C98", ls=(0, (4, 2)), lw=1.7),
    "grain": dict(color="#2E86C1", ls="-", lw=1.7),
    "rec": dict(color="#34495E", ls="-", lw=1.0),
    "chem": dict(color="#1B2631", ls="-", lw=1.5),
    "coll": dict(color="#E67E22", ls="-.", lw=1.5),
}
LEGEND = [("cr", "cosmic-ray ionization"),
          ("crphoto", "CR-induced photodissociation"),
          ("uv", "FUV photoreaction"),
          ("grain", "grain-assisted"),
          ("rec", "recombination (e)"),
          ("chem", "ion-molecule / neutral-neutral"),
          ("coll", "collisional (high T)")]

# (source, target, kind, label, arc[, label_t]). Label names the other
# reactant; label_t places it along the arrow (0 at source, 1 at target).
EDGES = [
    # cosmic rays
    ("H2", "H2+", "cr", "CR", 0.0),
    ("He", "He+", "cr", "CR", 0.25),
    ("H", "H+", "cr", "CR", 0.25),
    ("C", "C+", "cr", "CR", 0.25),
    ("CO", "C", "crphoto", "CR", 0.35),
    ("CO", "HCO+", "cr", "CR", 0.3),
    ("Si", "Si+", "cr", "CR", 0.3),
    # FUV
    ("C", "C+", "uv", "", -0.25),
    ("CHx", "C", "uv", "", 0.3),
    ("CO", "C", "uv", "", -0.25),
    ("OHx", "O", "uv", "", 0.3),
    ("H2", "H", "uv", "", 0.35),
    ("Si", "Si+", "uv", "", -0.3),
    # grains
    ("H", "H2", "grain", "H", 0.4),
    ("H+", "H", "grain", "e", -0.55),
    ("C+", "C", "grain", "e", -0.6),
    ("He+", "He", "grain", "e", 0.3),
    ("Si+", "Si", "grain", "e", 0.55),
    # recombination
    ("He+", "He", "rec", "e", 0.0),
    ("H3+", "H2", "rec", "e", -0.45),
    ("H3+", "H", "rec", "e", 0.2),
    ("C+", "C", "rec", "e", 0.0),
    ("HCO+", "CO", "rec", "e", 0.3),
    ("H+", "H", "rec", "e", 0.0),
    ("Si+", "Si", "rec", "e", 0.0),
    # ion-molecule and neutral-neutral
    ("H3+", "CHx", "chem", "C", 0.0),
    ("H3+", "OHx", "chem", "O", -0.25),
    ("CO", "HCO+", "chem", "H$_3^+$", -0.3),
    ("He+", "H+", "chem", "H$_2$", 0.12, 0.25),
    ("CO", "C+", "chem", "He$^+$", -0.3),
    ("C+", "CHx", "chem", "H$_2$", 0.0),
    ("C+", "HCO+", "chem", "OH$_x$", 0.0, 0.3),
    ("OHx", "HCO+", "chem", "C$^+$", 0.2, 0.72),
    ("CHx", "CO", "chem", "O", 0.25),
    ("OHx", "CO", "chem", "C", 0.25),
    ("H2+", "H3+", "chem", "H$_2$", 0.0),
    ("He+", "H2+", "chem", "H$_2$", 0.0),
    ("CHx", "C", "chem", "H", -0.3),
    ("OHx", "O", "chem", "O", -0.3),
    ("OHx", "O+", "chem", "He$^+$", 0.2),
    ("H2+", "H+", "chem", "H", 0.0),
    ("H+", "O+", "chem", "O", 0.12),
    ("O+", "H+", "chem", "H", 0.12),
    ("O+", "OHx", "chem", "H$_2$", 0.2),
    ("O+", "O", "chem", "H$_2$", 0.0),
    ("C+", "C", "chem", "H$_2$, e", 0.5),       # (22) C+ + H2 + e -> C + 2H
    ("H3+", "H2", "chem", "O, e", 0.45, 0.82),  # (24) H3+ + O + e -> H2 + O + H
    # collisional
    ("H2", "H", "coll", "H, H$_2$", 0.0),
    ("H", "H+", "coll", "e", -0.25),
]

# Two-body reactions drawn from one reactant to the product, whose other
# reactant would otherwise appear only as a text label. A dotted line joins that
# partner to the arrow, so He+ and H3+ show every reaction they drive.
PARTNERS = [("CO", "C+", "He+"), ("OHx", "O+", "He+"),
            ("CO", "HCO+", "H3+")]

BOX_W, BOX_H = 1.5, 0.8
LEGEND_FS, LEGEND_DY = 11.5, 0.56


def draw_arrow(ax, a, b, arc, st, shrink=2):
    """Arrow whose body carries the line pattern and whose head is solid.

    A dash or dot pattern applied to the whole patch breaks up the outline of
    the head too, so the body and the head are drawn as separate patches.
    """
    common = dict(connectionstyle=f"arc3,rad={arc}", color=st["color"],
                  shrinkA=shrink, shrinkB=shrink, zorder=1)
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", lw=st["lw"], ls=st["ls"],
                                 **common))
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>,head_length=6,head_width=3",
                                 lw=0, **common))


def box_edge(p, q):
    """Point where the segment p->q leaves the box centred on p."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    sx = (BOX_W / 2) / abs(dx) if dx else float("inf")
    sy = (BOX_H / 2) / abs(dy) if dy else float("inf")
    s = min(sx, sy)
    return (p[0] + s * dx, p[1] + s * dy)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", default="fig_gow17_network.png")
    args = parser.parse_args()

    # Limits hug the boxes and the legend graphics; legend text runs past the
    # right edge and bbox_inches="tight" keeps it.
    fig, ax = plt.subplots(figsize=(16, 7.3))
    ax.set_xlim(-2.7, 20.4)
    ax.set_ylim(-2.1, 9.3)
    ax.set_aspect("equal")
    ax.axis("off")

    # The coupled CO-HCO+ block.
    x0, y1 = POS["CO"][0], POS["CO"][1]
    ax.add_patch(FancyBboxPatch((x0 - 1.05, POS["HCO+"][1] - 0.65), 2.1,
                                y1 - POS["HCO+"][1] + 1.3,
                                boxstyle="round,pad=0.05,rounding_size=0.3",
                                fc="#FCF3CF", ec="#B7950B", lw=1.2, zorder=0))
    ax.text(x0 + 1.15, (y1 + POS["HCO+"][1]) / 2, "coupled\nblock", fontsize=9,
            color="#7D6608", va="center")

    def on_arc(a, b, arc, t):
        # arc3 is the quadratic Bezier with control point mid + rad*(dy, -dx).
        cx = (a[0] + b[0]) / 2 + arc * (b[1] - a[1])
        cy = (a[1] + b[1]) / 2 - arc * (b[0] - a[0])
        s = 1.0 - t
        return (s * s * a[0] + 2 * s * t * cx + t * t * b[0],
                s * s * a[1] + 2 * s * t * cy + t * t * b[1])

    midpoints = {}
    for src, dst, kind, label, arc, *rest in EDGES:
        label_t = rest[0] if rest else 0.5
        p, q = POS[src], POS[dst]
        a, b = box_edge(p, q), box_edge(q, p)
        midpoints[(src, dst, kind)] = on_arc(a, b, arc, 0.5)
        st = STYLE[kind]
        draw_arrow(ax, a, b, arc, st)
        # The arrow style already names CR, FUV, grain and recombination
        # processes; only the partner reactant of a two-body reaction is news.
        if label and kind in ("chem", "coll"):
            mx, my = on_arc(a, b, arc, label_t)
            ax.text(mx, my, label, fontsize=8, color=st["color"], ha="center",
                    va="center", zorder=3,
                    path_effects=[patheffects.withStroke(linewidth=4,
                                                         foreground="white")])

    for src, dst, partner in PARTNERS:
        m = midpoints[(src, dst, "chem")]
        start = box_edge(POS[partner], m)
        ax.plot([start[0], m[0]], [start[1], m[1]], ls=(0, (1, 2)), lw=1.3,
                color="#5D6D7E", zorder=1)
        ax.plot(*m, "o", ms=4, color="#5D6D7E", zorder=1.2)

    for name, (x, y) in POS.items():
        if name in NCR:
            ax.add_patch(FancyBboxPatch((x - BOX_W / 2 - 0.14, y - BOX_H / 2 - 0.14),
                                        BOX_W + 0.28, BOX_H + 0.28,
                                        boxstyle="round,pad=0.02,rounding_size=0.2",
                                        fc=NCR_COLOR[NCR[name]], ec="none", alpha=0.5,
                                        zorder=1.5))
        ghost = name in GHOSTS
        fc = "white" if ghost else ("#D6EAF8" if name in IONS else "#D5F5E3")
        ax.add_patch(FancyBboxPatch((x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
                                    boxstyle="round,pad=0.02,rounding_size=0.12",
                                    fc=fc, ec="#7F8C8D" if ghost else "#1B2631",
                                    ls="--" if ghost else "-",
                                    lw=2.8 if name in LIMITERS else 1.2, zorder=2))
        ax.text(x, y, TEX.get(name, name), ha="center", va="center", fontsize=13,
                style="italic" if ghost else "normal",
                color="#5D6D7E" if ghost else "black", zorder=4)
        if name in GS_ORDER:
            ax.text(x - BOX_W / 2 + 0.05, y + BOX_H / 2 - 0.02,
                    str(GS_ORDER.index(name) + 1), fontsize=7.5, ha="center",
                    va="center", color="white", zorder=5,
                    bbox=dict(boxstyle="circle,pad=0.18", fc="#1B2631", ec="none"))

    # Legend: arrow kinds, then box kinds.
    lx, ly = 19.3, 8.7
    for i, (kind, text) in enumerate(LEGEND):
        y = ly - LEGEND_DY * i
        st = STYLE[kind]
        draw_arrow(ax, (lx, y), (lx + 0.8, y), 0.0, st, shrink=0)
        ax.text(lx + 0.95, y, text, fontsize=LEGEND_FS, va="center")
    y = ly - LEGEND_DY * len(LEGEND)
    ax.plot([lx, lx + 0.8], [y, y], ls=(0, (1, 2)), lw=1.3, color="#5D6D7E")
    ax.plot(lx + 0.8, y, "o", ms=4, color="#5D6D7E")
    ax.text(lx + 0.95, y, "partner reactant (He$^+$, H$_3^+$)", fontsize=LEGEND_FS,
            va="center")
    y -= 0.7
    for fc, ec, ls, lw, text in [
            ("#D6EAF8", "#1B2631", "-", 1.2, "ion"),
            ("#D5F5E3", "#1B2631", "-", 1.2, "molecule"),
            ("white", "#7F8C8D", "--", 1.2, "ghost (from conservation laws)"),
            ("white", "#1B2631", "-", 2.8, "limits the substep")]:
        ax.add_patch(FancyBboxPatch((lx, y - 0.17), 0.8, 0.34,
                                    boxstyle="round,pad=0.01,rounding_size=0.06",
                                    fc=fc, ec=ec, ls=ls, lw=lw))
        ax.text(lx + 0.95, y, text, fontsize=LEGEND_FS, va="center")
        y -= LEGEND_DY
    for kind in ("ODE", "ghost", "equilibrium"):
        ax.add_patch(FancyBboxPatch((lx - 0.08, y - 0.25), 0.96, 0.5,
                                    boxstyle="round,pad=0.01,rounding_size=0.1",
                                    fc=NCR_COLOR[kind], ec="none", alpha=0.5))
        ax.text(lx + 0.95, y, NCR_LEGEND[kind], fontsize=LEGEND_FS, va="center")
        y -= LEGEND_DY
    ax.text(lx, y, "① order in block Gauss–Seidel", fontsize=LEGEND_FS,
            va="center")

    # Charge neutrality, as in GOW17Network::ComputeGhostSpecies_.
    ax.text(8.8, -1.0,
            r"$x_e = x({\rm He^+}) + x({\rm C^+}) + x({\rm HCO^+}) + x({\rm H_3^+})"
            r" + x({\rm H_2^+}) + x({\rm H^+}) + x({\rm O^+}) + x({\rm Si^+})$"
            "   (charge neutrality; electrons are a ghost species)",
            ha="center", va="center", fontsize=LEGEND_FS + 1)
    ax.text(8.8, -1.75,
            "Each of the 49 reactions is one arrow to its main product, except "
            "(15) H$_2$ + H and (16) H$_2$ + H$_2$, which share the collisional "
            "H$_2$ → H arrow, and C$^+$ + OH → HCO$^+$, drawn from both "
            "reactants.", ha="center", va="center", fontsize=LEGEND_FS - 1,
            color="#34495E")
    ax.set_title("GOW17 network: 12 integrated species, 5 ghost "
                 "species plus electrons", fontsize=13)
    out = pathlib.Path(args.out)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    print("wrote", out, "and", out.with_suffix(".pdf"))


if __name__ == "__main__":
    main()
