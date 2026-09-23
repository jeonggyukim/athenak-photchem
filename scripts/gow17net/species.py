"""The 18 species of the GOW17 network: 12 integrated, 6 ghosts."""
from dataclasses import dataclass


@dataclass(frozen=True)
class Species:
    name: str
    tex: str
    kind: str                 # "ion", "molecule" or "ghost"
    pseudo: bool = False      # CHx, OHx: lumped neutral and ionized hydrides
    gs_order: int | None = None      # position in the block Gauss-Seidel update
    limits_step: bool = False # limits the semi-implicit substep
    ncr: str | None = None           # how TIGRESS-NCR closes it, if it carries it


SPECIES = [
    Species("He+", "He$^+$", "ion", gs_order=1),
    Species("Si+", "Si$^+$", "ion", gs_order=2),
    Species("H2+", "H$_2^+$", "ion", gs_order=3),
    Species("H3+", "H$_3^+$", "ion", gs_order=4),
    Species("H+", "H$^+$", "ion", gs_order=5, ncr="non-equilibrium"),
    Species("O+", "O$^+$", "ion", gs_order=6, ncr="equilibrium"),
    Species("C+", "C$^+$", "ion", gs_order=7, ncr="equilibrium"),
    Species("CHx", "CH$_x$", "molecule", pseudo=True, gs_order=8),
    Species("OHx", "OH$_x$", "molecule", pseudo=True, gs_order=9),
    Species("CO", "CO", "molecule", gs_order=10, limits_step=True, ncr="equilibrium"),
    Species("HCO+", "HCO$^+$", "ion", gs_order=11),
    Species("H2", "H$_2$", "molecule", gs_order=12, limits_step=True,
            ncr="non-equilibrium"),
    Species("H", "H", "ghost", ncr="ghost"),
    Species("He", "He", "ghost"),
    Species("C", "C", "ghost", ncr="ghost"),
    Species("O", "O", "ghost", ncr="ghost"),
    Species("Si", "Si", "ghost"),
    Species("e-", "e$^-$", "ghost"),
]

BY_NAME = {s.name: s for s in SPECIES}
# Species drawn as boxes; electrons appear only as labels and in the x_e sum.
DRAWN = [s.name for s in SPECIES if s.name != "e-"]
