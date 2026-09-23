"""A plane-parallel PDR lit from the left: x is depth, log A_V, and each
species sits at the A_V where its equilibrium abundance peaks.

One row per element family; HCO+ goes in the O row, below CO, so the
CO-HCO+ cycle is a short vertical pair. The background is a cross-section of the slab:
zones named by what dominates, dashed lines at the H/H2, C+/C and C/CO fronts,
and the element fractions of H, H2, C+, C and CO in a strip above the rows.

The profile uses DESPOTIC's GOW17 at n_H = 1e3 cm^-3 and chi = 30. At chi = 1
the H/H2 front lies below A_V = 0.01, on top of the surface species; at
n_H = 1e2 or 30 cm^-3 with chi = 10 and xi_CR = 2e-16 s^-1, CO never exceeds
C (cosmic rays hold CO/C below 0.14), so there is no C/CO front. The profile
is cached in figures/gow17net/slab_profile.npz.
"""
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

from . import Layout
from .. import fluxes
from ..pdr_profile import OUTDIR, _crossing
from .. import arcs
from ..reactions import REACTIONS, index_of
from ..render import box_edge, on_arc

NH = 1.0e3
CHI = 30.0
T = 50.0
XI_CR = 2.0e-16
AV = np.logspace(-2, 1, 46)
CACHE = OUTDIR / "slab_profile.npz"
SPECIES = fluxes.gowmod.specListExtended


def _compute():
    x = np.empty((len(SPECIES), len(AV)))
    for j, av in enumerate(AV):
        net = fluxes.make_network(NH, T, CHI, av, XI_CR)
        xj = fluxes.equilibrium(net)
        x[:, j] = net.extendAbundances(xj, outdict=True).x.ravel()
    return x


def profile():
    """(A_V array, dict species -> equilibrium abundance per H nucleus)."""
    params = [NH, T, CHI, XI_CR]
    x = None
    if CACHE.exists():
        d = np.load(CACHE)
        if np.array_equal(d["AV"], AV) and list(d["species"]) == SPECIES \
                and np.array_equal(d["params"], params):
            x = d["x"]
    if x is None:
        x = _compute()
        OUTDIR.mkdir(parents=True, exist_ok=True)
        np.savez(CACHE, AV=AV, species=SPECIES, x=x, params=params)
    return AV.copy(), {s: x[i] for i, s in enumerate(SPECIES)}


def fronts(av, x):
    return {"H/H$_2$": _crossing(av, x["H"], 2 * x["H2"]),
            "C$^+$/C": _crossing(av, x["C+"], x["C"]),
            "C/CO": _crossing(av, x["C"], x["CO"])}


# Frame, in data units. x = X0 + K log10(A_V / A_V,min) across the slab.
X0, X1 = 0.0, 16.0
K = (X1 - X0) / np.log10(AV[-1] / AV[0])
ROWS = {"H": 8.2, "He": 6.75, "C": 5.3, "O": 3.85, "Si": 2.4}
FAMILY = {"H": "H", "H+": "H", "H2": "H", "H2+": "H", "H3+": "H",
          "He": "He", "He+": "He",
          "C": "C", "C+": "C", "CHx": "C", "CO": "C", "HCO+": "O",
          "O": "O", "O+": "O", "OHx": "O", "Si": "Si", "Si+": "Si"}
Y_ROW_LO, Y_ROW_HI = 1.6, 9.0       # zones span the rows
Y_STRIP = (9.3, 10.6)               # element-fraction curves
Y_AXIS = 1.35
BOX_W, BOX_H = 1.4, 0.7
GAP = 0.3           # least space between two boxes in a row
EDGE_GAP = 0.7      # space between boxes stacked at an edge of the slab


def av_to_x(av):
    return X0 + K * np.log10(np.asarray(av) / AV[0])


def _half_point(av, v, rising):
    """A_V where v last exceeds (falling) or first reaches (rising) half its max."""
    above = np.nonzero(v >= 0.5 * v.max())[0]
    return av[above[0]] if rising else av[above[-1]]


def _positions(av, x):
    """Species -> (x, y).

    A species whose abundance peaks inside the A_V grid sits at that A_V.
    One that peaks at either edge of the grid sits at that edge of the slab;
    several at one edge are ordered by where they fall to (or rise to) half
    their peak; on a tie, the one that keeps less of its peak at the far end
    of the slab sits nearer its own edge. A species within 10 per cent of its
    peak everywhere (He) takes the place after the surface species of the
    fullest row, still in the surface zone but clear of the arrows that
    leave those species. Boxes that overlap are then pushed apart
    within the slab.
    """
    rows = {}
    for s, fam in FAMILY.items():
        v = x[s]
        j = int(np.argmax(v))
        if v.min() >= 0.9 * v.max():
            key = ("flat", None)
        elif j == 0:
            key = ("surface", (_half_point(av, v, rising=False), v[-1] / v.max()))
        elif j == len(v) - 1:
            key = ("deep", (_half_point(av, v, rising=True), -v[0] / v.max()))
        else:
            key = ("peak", av[j])
        rows.setdefault(fam, []).append((s, key))
    step = BOX_W + GAP
    lo, hi = X0 + GAP + BOX_W / 2, X1 - GAP - BOX_W / 2
    n_surf = max(sum(k[0] == "surface" for _, k in items) for items in rows.values())
    pos = {}
    for fam, items in rows.items():
        xs = {}
        surf = sorted((k[1], s) for s, k in items if k[0] == "surface")
        deep = sorted(((k[1], s) for s, k in items if k[0] == "deep"), reverse=True)
        for i, (_, s) in enumerate(surf):
            xs[s] = lo + i * (BOX_W + EDGE_GAP)
        for i, (_, s) in enumerate(deep):
            xs[s] = hi - i * (BOX_W + EDGE_GAP)
        for s, k in items:
            if k[0] == "peak":
                xs[s] = float(av_to_x(k[1]))
            elif k[0] == "flat":
                xs[s] = lo + n_surf * (BOX_W + EDGE_GAP)
        order = sorted(xs, key=xs.get)
        for _ in range(200):
            moved = False
            for a, b in zip(order, order[1:]):
                d = step - (xs[b] - xs[a])
                if d > 1e-9:
                    xs[a] -= d / 2
                    xs[b] += d / 2
                    moved = True
            for s in order:
                xs[s] = min(max(xs[s], lo), hi)
            if not moved:
                break
        for s in order:
            pos[s] = (xs[s], ROWS[fam])
    return pos


def _curve(a, b, rad, n=60):
    t = np.linspace(0.0, 1.0, n)
    return np.array([on_arc(a, b, rad, ti) for ti in t])


def _arcs(pos, fixed=None):
    """rad per reaction, chosen so that curves miss the other boxes.

    Arrows between the same pair are handled together, in the order
    arcs.assign uses. Each takes the bulge (the curve's largest distance from
    the chord, in data units) that minimises boxes crossed, then leaving the
    row band, then its size; arrows of one pair keep bulges at least 0.3 apart.
    """
    fixed = fixed or {}
    groups = {}
    for i, r in enumerate(REACTIONS):
        pair, _ = arcs.canonical(r.src, r.dst)
        groups.setdefault(pair, []).append(i)
    rads = {}
    bulges = np.round(np.arange(-2.5, 2.51, 0.05), 2)
    pad_w, pad_h = BOX_W / 2 + 0.08, BOX_H / 2 + 0.08
    for (p, q), idx in groups.items():
        idx.sort(key=lambda i: (arcs._PROCESS_ORDER.get(REACTIONS[i].process, 99),
                                REACTIONS[i].athenak))
        a = box_edge(pos[p], pos[q], BOX_W, BOX_H)
        b = box_edge(pos[q], pos[p], BOX_W, BOX_H)
        length = np.hypot(b[0] - a[0], b[1] - a[1])
        others = [pos[s] for s in pos if s not in (p, q)]
        cost = {}
        for bu in bulges:
            rad = 2.0 * bu / length
            c = _curve(a, b, rad)
            hits = 0
            for ox, oy in others:
                inside = (np.abs(c[:, 0] - ox) < pad_w) & (np.abs(c[:, 1] - oy) < pad_h)
                hits += 30.0 * inside.any() + 0.5 * inside.sum()
            out = np.sum((c[:, 1] > Y_ROW_HI) | (c[:, 1] < Y_ROW_LO)
                         | (c[:, 0] < X0) | (c[:, 0] > X1 + 0.3))
            cost[bu] = hits + 3.0 * out + 1.5 * abs(bu)
        taken = []
        for i in idx:
            r = REACTIONS[i]
            _, forward = arcs.canonical(r.src, r.dst)
            if i in fixed:
                rads[i] = fixed[i]
                bu = (fixed[i] if forward else -fixed[i]) * length / 2
                taken.append(bu)
                continue
            ok = [bu for bu in bulges if all(abs(bu - t) >= 0.3 - 1e-9 for t in taken)]
            bu = min(ok, key=lambda v: (cost[v], abs(v)))
            taken.append(bu)
            rad = 2.0 * bu / length
            rads[i] = float(rad if forward else -rad)
    return rads


ZONE_COLORS = ["#F4ECF7", "#EAF2F8", "#FEF5E7", "#E8F6F3"]
ZONE_NAMES = ["H + C$^+$", "H$_2$ + C$^+$", "H$_2$ + C", "H$_2$ + CO"]
CURVES = [("H", 1, "H", "#7D3C98"), ("H2", 2, "H$_2$", "#2471A3"),
          ("C+", None, "C$^+$", "#A04000"), ("C", None, "C", "#B7950B"),
          ("CO", None, "CO", "#117A65")]
GREY = "#566573"


def _background(ax, layout):
    av, x = profile()
    fr = fronts(av, x)
    edges = [X0] + [float(av_to_x(a)) for a in fr.values()] + [X1]
    for (a, b), c, name in zip(zip(edges, edges[1:]), ZONE_COLORS, ZONE_NAMES):
        ax.add_patch(Rectangle((a, Y_ROW_LO), b - a, Y_STRIP[1] - Y_ROW_LO + 0.75,
                               fc=c, ec="none", zorder=0))
        ax.text((a + b) / 2, Y_STRIP[1] + 0.4, name, ha="center", va="center",
                fontsize=10.5, color="#1B2631", zorder=0.5)
    for (name, a), xf in zip(fr.items(), edges[1:-1]):
        ax.plot([xf, xf], [Y_AXIS, Y_STRIP[1] + 0.75], ls="--", lw=1.0,
                color=GREY, zorder=0.4)
        ax.text(xf, Y_STRIP[1] + 0.95, f"{name}\n$A_V$ = {a:.2f}", ha="center",
                va="bottom", fontsize=9, color=GREY, linespacing=1.1)

    # Element fractions of the main carriers, 0 to 1, in a strip above the rows.
    y0, y1 = Y_STRIP
    ax.plot([X0, X1], [y0, y0], color="#ABB2B9", lw=0.6, zorder=0.4)
    xc = x["C+"] + x["C"] + x["CO"] + x["CHx"] + x["HCO+"]
    xx = av_to_x(av)
    for k, (s, fac, label, c) in enumerate(CURVES):
        f = fac * x[s] if fac else x[s] / xc
        ax.plot(xx, y0 + (y1 - y0) * f, color=c, lw=1.6, alpha=0.8, zorder=0.6)
        ax.text(X1 + 0.12, y1 - 0.02 - 0.27 * k, label, color=c, fontsize=8.5,
                ha="left", va="center", zorder=0.7)
    ax.text(X0 - 0.15, (y0 + y1) / 2, "element\nfraction", ha="right",
            va="center", fontsize=8, color=GREY, linespacing=1.0)

    # Depth axis.
    ax.plot([X0, X1], [Y_AXIS, Y_AXIS], color=GREY, lw=1.0)
    for dec in (0.01, 0.1, 1, 10):
        for m in range(1, 10):
            a = dec * m
            if a < AV[0] * 0.999 or a > AV[-1] * 1.001:
                continue
            xt = float(av_to_x(a))
            major = m == 1
            ax.plot([xt, xt], [Y_AXIS, Y_AXIS - (0.16 if major else 0.08)],
                    color=GREY, lw=1.0 if major else 0.6)
            if major:
                ax.text(xt, Y_AXIS - 0.22, f"{a:g}", ha="center", va="top",
                        fontsize=9.5, color=GREY)
    ax.text((X0 + X1) / 2, Y_AXIS - 0.62,
            "$A_V$ [mag], depth from the illuminated surface", ha="center",
            va="top", fontsize=10, color=GREY)
    ax.text(X0, Y_AXIS - 0.62, "surface", ha="left", va="top", fontsize=9.5,
            color=GREY, style="italic")
    ax.text(X1, Y_AXIS - 0.62, "dark interior", ha="right", va="top",
            fontsize=9.5, color=GREY, style="italic")

    # FUV enters through the left face.
    for y in (ROWS["H"] - 0.75, (ROWS["C"] + ROWS["O"]) / 2, ROWS["Si"] - 0.4):
        ax.add_patch(FancyArrowPatch((X0 - 1.25, y), (X0 - 0.1, y),
                                     arrowstyle="-|>,head_length=7,head_width=4",
                                     lw=2.2, color="#7D3C98", zorder=0.5))
    ax.text(X0 - 0.7, (ROWS["C"] + ROWS["O"]) / 2 + 0.3, "FUV", ha="center",
            va="bottom", fontsize=11, color="#7D3C98", weight="bold")
    ax.plot([X0, X0], [Y_ROW_LO, Y_STRIP[1] + 0.75], color="#7D3C98", lw=1.5,
            zorder=0.5)
    for fam, y in ROWS.items():
        ax.text(X1 + 0.15, y, fam, ha="left", va="center", fontsize=10,
                color=GREY, style="italic")

    ax.text(X0 - 1.3, layout.ylim[1] - 0.1,
            "plane-parallel PDR, illuminated from one side:  "
            rf"$n_{{\rm H}} = 10^{{{np.log10(NH):.0f}}}$ cm$^{{-3}}$, "
            rf"$\chi = {CHI:g}$, $\xi_{{\rm CR}} = 2\times10^{{-16}}$ s$^{{-1}}$, "
            rf"$T = {T:.0f}$ K, equilibrium", ha="left", va="top", fontsize=11,
            color="#1B2631")


_av, _x = profile()
_POS = _positions(_av, _x)

# CO -> C+ (partner He+) bends down, so the dotted line from He+ to its
# midpoint does not run along the He row through the He box.
# H3+ + O, e -> H2 runs straight, so its label sits in the H2-H3+ gap.
_FIXED = {index_of("2body", 4): -0.16, index_of("2body", 24): 0.0}
_ARCS = _arcs(_POS, _FIXED)
# He+ + H2 -> H+ and H + e -> H+ would put their labels on top of each other.
_LABELS = {index_of("2body", 3): 0.4, index_of("2body", 17): 0.75,
           index_of("2body", 6, "OHx"): 0.35}

SLAB = Layout("slab", _POS, xlim=(-1.6, 23.6), ylim=(-0.9, 12.9),
              legend_xy=(17.4, 10.6), equation_xy=(X0 + (X1 - X0) / 2, -0.15),
              caption_xy=(X0 + (X1 - X0) / 2, -0.7), box_w=BOX_W, box_h=BOX_H,
              background=_background, arc_overrides=_ARCS,
              label_overrides=_LABELS)
