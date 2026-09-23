"""The grid layout of the committed network figure, and a variant with rows
grouped by charge (neutral row on top, ion row below)."""
from . import Layout
from ..reactions import index_of

_POS = {
    "Si": (-1.8, 8.0), "Si+": (-1.8, 5.5),
    "H": (1.0, 8.0), "H2": (4.2, 8.0), "H+": (1.0, 3.0),
    "H2+": (4.2, 5.5), "H3+": (4.2, 3.0),
    "He": (7.4, 8.0), "He+": (7.4, 5.5),
    "C": (10.6, 8.0), "C+": (13.8, 8.0), "CO": (17.0, 8.0),
    "CHx": (11.0, 5.5), "HCO+": (17.0, 5.0),
    "OHx": (13.8, 3.0), "O+": (10.6, 0.6), "O": (13.8, 0.6),
}

_FRAME = dict(xlim=(-2.9, 26.4), ylim=(-2.5, 9.4), legend_xy=(19.4, 8.9),
              equation_xy=(8.6, -1.2), caption_xy=(8.6, -2.0))

# The two H3+ -> H2 arrows pass either side of H2+, which sits between them.
_ARCS = {index_of("2body", 10): -0.45, index_of("2body", 24): 0.45,
         # H+ <-> O+ charge exchange: the straight line already passes below
         # H3+, so O+ -> H+ stays nearly straight and H+ -> O+ bends under it.
         index_of("2body", 28): -0.05, index_of("2body", 27): 0.3}
# He+ + H2 -> H+ crosses H2+ -> H3+ where both labels would sit.
_LABELS = {index_of("2body", 24): 0.82, index_of("2body", 3): 0.25,
           index_of("2body", 6, "C+"): 0.3, index_of("2body", 6, "OHx"): 0.72}

GRID = Layout("grid", dict(_POS), arc_overrides=_ARCS, label_overrides=_LABELS,
              **_FRAME)

# C+ and CHx swap places: the top row then holds only neutrals and the second
# row only ions, apart from the oxygen block.
_ROWS = dict(_POS)
_ROWS["C+"], _ROWS["CHx"] = (11.0, 5.5), (13.8, 8.0)
GRID_CHARGE_ROWS = Layout("grid_charge_rows", _ROWS, arc_overrides=_ARCS,
                          label_overrides=_LABELS, **_FRAME)
