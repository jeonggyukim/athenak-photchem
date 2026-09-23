"""Layouts: where each species box sits, and the fixed frame around them.

A layout is data. Every figure drawn with one layout uses the same figure
size, axes rectangle and data limits, so boxes sit at the same pixel position
in every variant and slides can be flipped without the diagram moving.
"""
from dataclasses import dataclass, field


@dataclass
class Layout:
    name: str
    positions: dict                 # species -> (x, y) in data units
    xlim: tuple
    ylim: tuple
    legend_xy: tuple                # top-left of the legend, data units
    equation_xy: tuple              # centre of the x_e equation
    caption_xy: tuple               # centre of the caption line
    box_w: float = 1.5
    box_h: float = 0.8
    figsize: tuple = (13.333, 7.5)  # a 16:9 slide
    dpi: int = 200
    axes_rect: tuple = (0.01, 0.02, 0.98, 0.90)  # figure fraction, below the title
    background: object = None       # callable(ax, layout) or None
    arc_overrides: dict = field(default_factory=dict)    # reaction index -> rad
    label_overrides: dict = field(default_factory=dict)  # reaction index -> t


def get(name):
    from . import grid
    layouts = {"grid": grid.GRID, "grid_charge_rows": grid.GRID_CHARGE_ROWS}
    return layouts[name]
