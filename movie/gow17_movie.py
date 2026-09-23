"""ManimGL movie: the GOW17 network built up one physical process at a time.

Species, reactions, build steps and captions all come from gow17net, so the
movie cannot drift from the figures. The layout is chosen by the environment
variable GOW17_LAYOUT (default "grid").

Render (conda env manimgl, from the repository root):
    manimgl movie/gow17_movie.py BuildNetwork -w -l     # 480p preview
    manimgl movie/gow17_movie.py BuildNetwork -w --hd   # 1080p
"""
import math

import numpy as np
import os
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

from manimlib import (DOWN, LEFT, RIGHT, UP, WHITE, AnimationGroup,  # noqa: E402
                      ArcBetweenPoints, ArrowTip, DashedVMobject, Dot, FadeIn,
                      FadeOut, FadeTransform, Indicate, LaggedStart, MoveAlongPath,
                      RoundedRectangle, Scene, ShowCreation, Succession, Tex,
                      TexText, VGroup)

from gow17net import arcs, layouts, style  # noqa: E402
from gow17net.builds import CAPTIONS  # noqa: E402
from gow17net.reactions import REACTIONS, STEPS  # noqa: E402
from gow17net.render import box_edge  # noqa: E402
from gow17net.species import BY_NAME  # noqa: E402

LAYOUT = layouts.get(os.environ.get("GOW17_LAYOUT", "grid"))
BG = "#101418"
PAST = "#4D5656"
PROCESS_COLOR = dict((k, v["color"]) for k, v in style.PROCESS_STYLE.items())
PROCESS_COLOR.update(chem="#F4F6F7", rec="#AEB6BF")      # visible on dark
LONG_ARROW = 5.0          # data units; arrows longer than this bend less
ZONE_TINT = ["#8E44AD", "#2E86C1", "#D68910", "#17A589"]  # H+C+, H2+C+, H2+C, H2+CO
FILL = {"ion": "#D6EAF8", "molecule": "#D5F5E3", "ghost": BG}


def _frame_map():
    """Map layout data coordinates onto the manim frame (14.2 x 8 units),
    keeping the top band for the title and the bottom band for captions."""
    xs = [p[0] for p in LAYOUT.positions.values()]
    ys = [p[1] for p in LAYOUT.positions.values()]
    x0, x1 = min(xs) - LAYOUT.box_w, max(xs) + LAYOUT.box_w
    y0, y1 = min(ys) - LAYOUT.box_h, max(ys) + LAYOUT.box_h
    band = 4.6 if LAYOUT.name == "slab" else 5.4
    s = min(13.0 / (x1 - x0), band / (y1 - y0))
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2

    def to_frame(p):
        return ((p[0] - cx) * s, (p[1] - cy) * s + 0.35, 0.0)
    return to_frame, s


TO_FRAME, SCALE = _frame_map()


DASHED = ("uv", "crphoto", "coll")


def _make_arrow(start, end, angle, color, process):
    """An arrow as a stroked arc (the body) plus a separate tip.

    ShowCreation on ManimGL's filled Arrow traces its outline, which shows up
    as a thin wedge with the tip first; a stroked body draws cleanly from the
    source to the target, and the tip is added when the body arrives.
    """
    start, end = np.array(start, dtype=float), np.array(end, dtype=float)
    body = ArcBetweenPoints(start, end, angle=angle if abs(angle) > 1e-6 else 1e-6)
    body.set_stroke(color, width=3.0)
    path = body.copy()
    if process in DASHED:
        body = DashedVMobject(body, num_dashes=max(6, int(body.get_arc_length() / 0.12)))
        body.set_stroke(color, width=3.0)
    tangent = path.get_end() - path.point_from_proportion(0.96)
    tip = ArrowTip(angle=math.atan2(tangent[1], tangent[0]), width=0.12, length=0.14)
    tip.set_fill(color, opacity=1.0).move_to(path.get_end() - 0.05 * tangent /
                                             (1e-9 + math.hypot(tangent[0], tangent[1])))
    return VGroup(body, tip), path


def _slab_background():
    """Zones, fronts, A_V axis and FUV arrows for the plane-parallel PDR."""
    from manimlib import Arrow as _Arrow, DashedLine, Rectangle, Text

    from gow17net.layouts import slab
    av, x = slab.profile()
    fronts = slab.fronts(av, x)
    edges = [slab.X0] + [float(slab.av_to_x(a)) for a in fronts.values()] + [slab.X1]
    y0, y1 = slab.Y_ROW_LO, slab.Y_ROW_HI
    zones, labels, lines = VGroup(), VGroup(), VGroup()
    for xa, xb, color, name in zip(edges[:-1], edges[1:], ZONE_TINT, slab.ZONE_NAMES):
        p0, p1 = TO_FRAME((xa, y0)), TO_FRAME((xb, y1))
        rect = Rectangle(width=p1[0] - p0[0], height=p1[1] - p0[1])
        rect.set_fill(color, opacity=0.22).set_stroke(width=0)
        rect.move_to(((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2, 0))
        zones.add(rect)
        lab = Tex(r"\mathrm{" + name.replace("$", "") + "}").set_color("#AEB6BF")
        lab.set_height(0.18).move_to(rect.get_top() + 0.16 * DOWN)
        labels.add(lab)
    for (name, a), xf in zip(fronts.items(), edges[1:-1]):
        top, bot = TO_FRAME((xf, y1)), TO_FRAME((xf, y0))
        lines.add(DashedLine(bot, top, dash_length=0.08).set_stroke("#AEB6BF", 1.5))
        tag = Tex(r"A_V=" + f"{a:.2f}").set_color("#AEB6BF").set_height(0.16)
        lines.add(tag.next_to(DashedLine(bot, top), DOWN, buff=0.06))
    fuv = VGroup()
    for yy in (y0 + 1.2, (y0 + y1) / 2, y1 - 1.2):
        p = TO_FRAME((slab.X0, yy))
        fuv.add(_Arrow((p[0] - 0.9, p[1], 0), (p[0] - 0.05, p[1], 0), buff=0,
                       fill_color="#AF7AC5", thickness=3))
    fuv_label = Text("FUV", font_size=24).set_color("#AF7AC5")
    fuv_label.next_to(fuv[-1], UP, buff=0.1)
    axis_l = Text("surface", font_size=18).set_color("#AEB6BF")
    axis_r = Text("dark interior", font_size=18).set_color("#AEB6BF")
    axis_l.next_to(TO_FRAME((slab.X0, y0)), DOWN, buff=0.3).shift(0.4 * RIGHT)
    axis_r.next_to(TO_FRAME((slab.X1, y0)), DOWN, buff=0.3).shift(0.6 * LEFT)
    return VGroup(zones, labels, lines), VGroup(fuv, fuv_label, axis_l, axis_r)


class BuildNetwork(Scene):
    def setup(self):
        self.camera.background_rgba = [0.063, 0.078, 0.094, 1.0]

    def construct(self):
        w, h = LAYOUT.box_w * SCALE, LAYOUT.box_h * SCALE
        boxes = {}
        for name, p in LAYOUT.positions.items():
            sp = BY_NAME[name]
            box = RoundedRectangle(width=w, height=h, corner_radius=0.08)
            box.set_fill(FILL[sp.kind], opacity=1.0)
            box.set_stroke("#7F8C8D" if sp.kind == "ghost" else "#1B2631", 2)
            box.move_to(TO_FRAME(p))
            label = Tex(r"\mathrm{" + sp.tex.replace("$", "") + "}"
                        if "$" in sp.tex else r"\mathrm{" + sp.tex + "}")
            label.set_color("#AAB7B8" if sp.kind == "ghost" else "#000000")
            label.set_height(min(h * 0.55, 0.32)).move_to(box)
            boxes[name] = VGroup(box, label)

        rads = arcs.assign(REACTIONS, overrides=LAYOUT.arc_overrides)
        arrows = {}
        labels = {}
        for i, (r, rad) in enumerate(zip(REACTIONS, rads)):
            a = box_edge(LAYOUT.positions[r.src], LAYOUT.positions[r.dst],
                         LAYOUT.box_w, LAYOUT.box_h)
            b = box_edge(LAYOUT.positions[r.dst], LAYOUT.positions[r.src],
                         LAYOUT.box_w, LAYOUT.box_h)
            # matplotlib arc3 with rad r has sagitta r * chord / 2; a circular
            # arc through the same points has path_arc = 4 atan(r).
            # The bulge of an arc grows with its length; scale the bend of long
            # arrows down so they stay inside the frame (same factor for all
            # arrows of one pair, so their order is kept).
            chord = math.dist(a, b)
            rad_eff = (rad * min(1.0, LONG_ARROW / chord)
                       if chord > 0 and LAYOUT.name == "slab" else rad)
            arrows[i] = _make_arrow(TO_FRAME(a), TO_FRAME(b), -4.0 * math.atan(rad_eff),
                                    PROCESS_COLOR[r.process], r.process)
            # Partner reactant of a two-body reaction, at the label position
            # the layout chose (fraction t along the drawn arc).
            if r.partner and r.process in ("chem", "coll", "ncr"):
                t = LAYOUT.label_overrides.get(i, 0.5)
                lab = Tex(r"\mathrm{" + r.partner.replace("$", "") + "}")
                lab.set_color(PROCESS_COLOR[r.process]).set_height(0.14)
                lab.set_backstroke(BG, 5)
                lab.move_to(arrows[i][1].point_from_proportion(t))
                labels[i] = lab

        if LAYOUT.name == "slab":
            zones, light = _slab_background()
            intro = TexText(r"A plane-parallel PDR, lit from the left").set_color(WHITE)
            intro.scale(0.8).to_edge(UP)
            self.play(FadeIn(intro), FadeIn(zones), run_time=2)
            self.play(LaggedStart(*[ShowCreation(a) for a in light[0]], lag_ratio=0.3),
                      FadeIn(light[1:]), run_time=2)
            self.wait(2)
            self.play(FadeOut(intro))

        title = TexText(STEPS[1]).set_color(WHITE).scale(0.8).to_edge(UP)
        caption = self._caption(1)
        self.play(FadeIn(title), LaggedStart(*[FadeIn(b) for b in boxes.values()],
                                             lag_ratio=0.05), run_time=3)
        self.play(FadeIn(caption))
        self.wait(2)

        shown = []
        for step in range(2, max(STEPS) + 1):
            new = [i for i, r in enumerate(REACTIONS) if r.step == step]
            if step == 11 or not new:
                continue
            new_title = TexText(f"{step}. {STEPS[step]}").set_color(WHITE)
            new_title.scale(0.8).to_edge(UP)
            new_caption = self._caption(step)
            fades = []
            for i in shown:
                body, tip = arrows[i][0]
                fades += [body.animate.set_stroke(PAST), tip.animate.set_fill(PAST)]
                if i in labels:
                    fades.append(labels[i].animate.set_fill(PAST))
            self.play(FadeTransform(title, new_title), FadeOut(caption),
                      *fades, run_time=1.0)
            beats = []
            for i in new:
                (body, tip), path = arrows[i]
                color = PROCESS_COLOR[REACTIONS[i].process]
                dot = Dot(radius=0.06).set_fill(color).move_to(path.get_start())
                extra = []
                label = labels.get(i)
                if label is not None:
                    extra.append(FadeIn(label, run_time=0.3))
                beats.append(Succession(
                    ShowCreation(body, run_time=0.8),
                    AnimationGroup(FadeIn(tip, run_time=0.2), *extra,
                                   MoveAlongPath(dot, path, run_time=0.6),
                                   Indicate(boxes[REACTIONS[i].dst], scale_factor=1.12,
                                            color=color, run_time=0.6)),
                    FadeOut(dot, run_time=0.1)))
            self.play(LaggedStart(*beats, lag_ratio=0.35),
                      run_time=max(2.5, 0.55 * len(new) + 1.2))
            self.play(FadeIn(new_caption))
            self.wait(3)
            title, caption = new_title, new_caption
            shown += new
        self.wait(2)

    def _caption(self, step):
        text = CAPTIONS[step].replace("–", "--")
        cap = TexText(text).set_color(WHITE)
        cap.set_width(min(cap.get_width(), 13.0)).to_edge(DOWN, buff=0.35)
        return cap
