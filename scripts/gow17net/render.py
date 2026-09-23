"""Draw the network for one layout and one variant.

The frame is fixed per layout (figure size, dpi, axes rectangle, data limits,
equal aspect) and figures are saved without cropping, so every variant of a
layout places every box at the same pixel position.

Variants:
- `step`: build-sequence step (reactions.STEPS). Reactions introduced at this
  step are drawn in colour, earlier ones grey, later ones fully transparent so
  they keep their place. None draws every GOW17 reaction in colour.
- `widths`: {reaction index: line width}, e.g. from widths.py.
- `annotations`: Gauss-Seidel numbers, step-limiter borders and NCR halos.
- `theme`: "light" or "dark".
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib import patheffects  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch  # noqa: E402

from . import arcs, style  # noqa: E402
from .reactions import REACTIONS, STEPS  # noqa: E402
from .species import BY_NAME  # noqa: E402

THEMES = {
    "light": dict(bg="white", fg="black", ghost_text="#5D6D7E"),
    "dark": dict(bg="#101418", fg="#ECF0F1", ghost_text="#AAB7B8"),
}
# Partner reactants that are species boxes get a dotted line to the arrow.
PARTNER_SPECIES = {"He$^+$": "He+", "H$_3^+$": "H3+"}
ANNOTATION_STEP = 11
NCR_STEP = 12


def box_edge(p, q, w, h):
    """Point where the segment p->q leaves the w x h box centred on p."""
    dx, dy = q[0] - p[0], q[1] - p[1]
    sx = (w / 2) / abs(dx) if dx else float("inf")
    sy = (h / 2) / abs(dy) if dy else float("inf")
    s = min(sx, sy)
    return (p[0] + s * dx, p[1] + s * dy)


def on_arc(a, b, rad, t):
    """Point at parameter t on matplotlib's arc3 curve from a to b."""
    cx = (a[0] + b[0]) / 2 + rad * (b[1] - a[1])
    cy = (a[1] + b[1]) / 2 - rad * (b[0] - a[0])
    s = 1.0 - t
    return (s * s * a[0] + 2 * s * t * cx + t * t * b[0],
            s * s * a[1] + 2 * s * t * cy + t * t * b[1])


def _arrow(ax, a, b, rad, color, lw, ls, alpha, zorder=1):
    # Body and head are separate patches, so a dash pattern does not break up
    # the outline of the head.
    common = dict(connectionstyle=f"arc3,rad={rad}", color=color, alpha=alpha,
                  shrinkA=2, shrinkB=2, zorder=zorder)
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-", lw=lw, ls=ls, **common))
    ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>,head_length=6,head_width=3",
                                 lw=0, **common))


def _reaction_state(r, step):
    """'current', 'past' or 'hidden' for a reaction at a build step."""
    if step is None:
        return "current" if r.network == "gow17" else "hidden"
    if r.step == step:
        return "current"
    return "past" if r.step < step else "hidden"


def draw(layout, step=None, widths=None, annotations=None, theme="light",
         title=None, caption=None, legend=True):
    """Return (figure, {species: box patch})."""
    th = THEMES[theme]
    if annotations is None:
        annotations = step is None or step >= ANNOTATION_STEP
    fig = plt.figure(figsize=layout.figsize, dpi=layout.dpi, facecolor=th["bg"])
    ax = fig.add_axes(layout.axes_rect)
    ax.set_xlim(*layout.xlim)
    ax.set_ylim(*layout.ylim)
    ax.set_aspect("equal", adjustable="box", anchor="C")
    ax.axis("off")
    ax.set_facecolor(th["bg"])
    if layout.background is not None:
        layout.background(ax, layout)

    pos, w, h = layout.positions, layout.box_w, layout.box_h
    rads = arcs.assign(REACTIONS, overrides=layout.arc_overrides)
    stroke = [patheffects.withStroke(linewidth=4, foreground=th["bg"])]

    for i, (r, rad) in enumerate(zip(REACTIONS, rads)):
        state = _reaction_state(r, step)
        if state == "hidden":
            continue
        a = box_edge(pos[r.src], pos[r.dst], w, h)
        b = box_edge(pos[r.dst], pos[r.src], w, h)
        st = style.PROCESS_STYLE[r.process]
        color = st["color"] if state == "current" else style.PAST_COLOR
        lw = widths.get(i, st["lw"]) if widths else st["lw"]
        _arrow(ax, a, b, rad, color, lw, st["ls"], 1.0)
        if r.partner and r.process in ("chem", "coll", "ncr"):
            t = layout.label_overrides.get(i, 0.5)
            mx, my = on_arc(a, b, rad, t)
            ax.text(mx, my, r.partner, fontsize=style.FONT["label"], color=color,
                    ha="center", va="center", zorder=3, path_effects=stroke)
            partner = PARTNER_SPECIES.get(r.partner)
            if partner and partner != r.src:
                m = on_arc(a, b, rad, 0.5)
                s = box_edge(pos[partner], m, w, h)
                ax.plot([s[0], m[0]], [s[1], m[1]], ls=(0, (1, 2)), lw=1.3,
                        color=style.PARTNER_COLOR if state == "current"
                        else style.PAST_COLOR, zorder=1)
                ax.plot(*m, "o", ms=4, color=style.PARTNER_COLOR, zorder=1.2)

    boxes = {}
    for name, (x, y) in pos.items():
        sp = BY_NAME[name]
        if annotations and sp.ncr:
            ax.add_patch(FancyBboxPatch((x - w / 2 - 0.14, y - h / 2 - 0.14),
                                        w + 0.28, h + 0.28,
                                        boxstyle="round,pad=0.02,rounding_size=0.2",
                                        fc=style.NCR_HALO[sp.ncr], ec="none",
                                        alpha=0.5, zorder=1.5))
        boxes[name] = ax.add_patch(FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.02,rounding_size=0.12",
            fc=style.BOX_FILL[sp.kind], ec=style.BOX_EDGE[sp.kind],
            ls="--" if sp.kind == "ghost" else "-",
            lw=2.8 if (annotations and sp.limits_step) else 1.2, zorder=2))
        ax.text(x, y, sp.tex, ha="center", va="center", fontsize=style.FONT["species"],
                style="italic" if sp.kind == "ghost" else "normal",
                color=th["ghost_text"] if sp.kind == "ghost" else "black", zorder=4)
        if annotations and sp.gs_order:
            ax.text(x - w / 2 + 0.05, y + h / 2 - 0.02, str(sp.gs_order),
                    fontsize=7.5, ha="center", va="center", color="white", zorder=5,
                    bbox=dict(boxstyle="circle,pad=0.18", fc="#1B2631", ec="none"))

    if legend:
        _legend(ax, layout, th, annotations, step)
    ax.text(*layout.equation_xy,
            r"$x_e = x({\rm He^+}) + x({\rm C^+}) + x({\rm HCO^+}) + x({\rm H_3^+})"
            r" + x({\rm H_2^+}) + x({\rm H^+}) + x({\rm O^+}) + x({\rm Si^+})$",
            ha="center", va="center", fontsize=style.FONT["legend"], color=th["fg"])
    if caption:
        ax.text(*layout.caption_xy, caption, ha="center", va="center",
                fontsize=style.FONT["caption"], color=th["fg"])
    if title is None:
        title = "GOW17 network" + (f": {STEPS[step]}" if step else "")
    fig.suptitle(title, fontsize=style.FONT["title"], color=th["fg"], y=0.97)
    return fig, boxes


def _legend(ax, layout, th, annotations, step):
    lx, ly = layout.legend_xy
    dy = 0.55
    fs = style.FONT["legend"] - 1
    kinds = ["cr", "crphoto", "uv", "grain", "rec", "chem", "coll"]
    if step is not None and step >= NCR_STEP:
        kinds.append("ncr")
    y = ly
    for k in kinds:
        st = style.PROCESS_STYLE[k]
        _arrow(ax, (lx, y), (lx + 0.8, y), 0.0, st["color"], st["lw"], st["ls"], 1.0)
        ax.text(lx + 0.95, y, style.PROCESS_LABEL[k], fontsize=fs, va="center",
                color=th["fg"])
        y -= dy
    y -= 0.15
    for kind, text in (("ion", "ion"), ("molecule", "molecule"),
                       ("ghost", "ghost (from conservation laws)")):
        ax.add_patch(FancyBboxPatch((lx, y - 0.17), 0.8, 0.34,
                                    boxstyle="round,pad=0.01,rounding_size=0.06",
                                    fc=style.BOX_FILL[kind], ec=style.BOX_EDGE[kind],
                                    ls="--" if kind == "ghost" else "-", lw=1.2))
        ax.text(lx + 0.95, y, text, fontsize=fs, va="center", color=th["fg"])
        y -= dy
    if annotations:
        ax.add_patch(FancyBboxPatch((lx, y - 0.17), 0.8, 0.34,
                                    boxstyle="round,pad=0.01,rounding_size=0.06",
                                    fc="white", ec="#1B2631", lw=2.8))
        ax.text(lx + 0.95, y, "limits the substep", fontsize=fs, va="center",
                color=th["fg"])
        y -= dy
        for key in ("non-equilibrium", "ghost", "equilibrium"):
            ax.add_patch(FancyBboxPatch((lx - 0.08, y - 0.25), 0.96, 0.5,
                                        boxstyle="round,pad=0.01,rounding_size=0.1",
                                        fc=style.NCR_HALO[key], ec="none", alpha=0.5))
            ax.text(lx + 0.95, y, style.NCR_LABEL[key], fontsize=fs, va="center",
                    color=th["fg"])
            y -= dy
        ax.text(lx, y, "① order in block Gauss–Seidel", fontsize=fs,
                va="center", color=th["fg"])


def save(fig, out):
    """Save without cropping, so the frame is identical across variants."""
    fig.savefig(out, dpi=fig.dpi, facecolor=fig.get_facecolor())
    plt.close(fig)
