"""The 49 reactions of AthenaK's GOW17 network, one entry per drawn arrow.

`athenak` holds the reaction ids as numbered within each group of
GOW17Network::CDRates in src/chemistry/network/gow17.hpp. An arrow runs from one
reactant (`src`) to the main product (`dst`); `partner` labels the other
reactant. `despotic` holds the matching keys in DESPOTIC's GOW17 tables, used
for arrow widths. `step` is the slide build step that introduces the arrow.

Two exceptions to one arrow per reaction: 2-body (15) H2 + H and (16) H2 + H2
share the collisional H2 -> H arrow, and 2-body (6) C+ + OH -> HCO+ is drawn
from both reactants (the OHx arrow is a second entry with the same id).
"""
from dataclasses import dataclass

from .species import BY_NAME

PROCESSES = ("cr", "crphoto", "uv", "grain", "rec", "chem", "coll")

# Build steps, in slide order.
STEPS = {
    1: "species",
    2: "cosmic-ray ionization",
    3: "FUV photoreactions",
    4: "grain processes",
    5: "radiative recombination",
    6: "hydrogen ion chemistry",
    7: "carbon chemistry",
    8: "oxygen chemistry",
    9: "CO-HCO+ cycle",
    10: "collisional (high T)",
    11: "solver annotations",
    12: "processes only in NCR",
}


@dataclass(frozen=True)
class Reaction:
    athenak: tuple            # ((group, number), ...)
    src: str
    dst: str
    process: str
    step: int
    partner: str | None = None
    despotic: tuple = ()      # (("cr"|"ph"|"tb", index), ...)
    network: str = "gow17"    # "gow17", or "ncr" for NCR-only processes


def _r(ids, src, dst, process, step, partner=None, despotic=()):
    ids = ids if isinstance(ids[0], tuple) else (ids,)
    return Reaction(ids, src, dst, process, step, partner, tuple(despotic))


REACTIONS = [
    # cosmic rays
    _r(("cr", 0), "H2", "H2+", "cr", 2, despotic=[("cr", "H2")]),
    _r(("cr", 1), "He", "He+", "cr", 2, despotic=[("cr", 0)]),
    _r(("cr", 2), "H", "H+", "cr", 2, despotic=[("cr", "H")]),
    _r(("cr", 3), "C", "C+", "cr", 2, despotic=[("cr", 1), ("cr", 3)]),
    _r(("cr", 4), "CO", "C", "crphoto", 3, despotic=[("cr", 4)]),
    _r(("cr", 5), "CO", "HCO+", "cr", 9, despotic=[("cr", 2)]),
    _r(("cr", 6), "Si", "Si+", "cr", 2, despotic=[("cr", 5)]),
    # FUV
    _r(("photo", 0), "C", "C+", "uv", 3, despotic=[("ph", 0)]),
    _r(("photo", 1), "CHx", "C", "uv", 3, despotic=[("ph", 1)]),
    _r(("photo", 2), "CO", "C", "uv", 3, despotic=[("ph", 2)]),
    _r(("photo", 3), "OHx", "O", "uv", 3, despotic=[("ph", 3)]),
    _r(("photo", 4), "H2", "H", "uv", 3, despotic=[("ph", 5)]),
    _r(("photo", 5), "Si", "Si+", "uv", 3, despotic=[("ph", 4)]),
    # grains
    _r(("grain", 0), "H", "H2", "grain", 4, despotic=[("tb", 31)]),
    _r(("grain", 1), "H+", "H", "grain", 4, despotic=[("tb", 32)]),
    _r(("grain", 2), "C+", "C", "grain", 4, despotic=[("tb", 33)]),
    _r(("grain", 3), "He+", "He", "grain", 4, despotic=[("tb", 34)]),
    _r(("grain", 4), "Si+", "Si", "grain", 4, despotic=[("tb", 35)]),
    # radiative recombination
    _r(("2body", 9), "He+", "He", "rec", 5, despotic=[("tb", 13)]),
    _r(("2body", 10), "H3+", "H2", "rec", 5, despotic=[("tb", 14)]),
    _r(("2body", 11), "C+", "C", "rec", 5, despotic=[("tb", 16)]),
    _r(("2body", 14), "H+", "H", "rec", 5, despotic=[("tb", 20)]),
    _r(("2body", 18), "H3+", "H", "rec", 5, despotic=[("tb", 15)]),
    _r(("2body", 23), "Si+", "Si", "rec", 5, despotic=[("tb", 27)]),
    # hydrogen ion chemistry
    _r(("2body", 3), "He+", "H+", "chem", 6, "H$_2$", [("tb", 6)]),
    _r(("2body", 13), "H2+", "H3+", "chem", 6, "H$_2$", [("tb", 18)]),
    _r(("2body", 19), "He+", "H2+", "chem", 6, "H$_2$", [("tb", 24)]),
    _r(("2body", 26), "H2+", "H+", "chem", 6, "H", [("tb", 19)]),
    # carbon chemistry
    _r(("2body", 0), "H3+", "CHx", "chem", 7, "C", [("tb", 0)]),
    _r(("2body", 4), "CO", "C+", "chem", 7, "He$^+$", [("tb", 7)]),
    _r(("2body", 5), "C+", "CHx", "chem", 7, "H$_2$", [("tb", 8)]),
    _r(("2body", 7), "CHx", "CO", "chem", 7, "O", [("tb", 11)]),
    _r(("2body", 20), "CHx", "C", "chem", 7, "H", [("tb", 25)]),
    _r(("2body", 22), "C+", "C", "chem", 7, "H$_2$, e", [("tb", 9)]),
    # oxygen chemistry
    _r(("2body", 1), "H3+", "OHx", "chem", 8, "O", [("tb", 1)]),
    _r(("2body", 8), "OHx", "CO", "chem", 8, "C", [("tb", 12)]),
    _r(("2body", 21), "OHx", "O", "chem", 8, "O", [("tb", 26)]),
    _r(("2body", 24), "H3+", "H2", "chem", 8, "O, e", [("tb", 2)]),
    _r(("2body", 25), "OHx", "O+", "chem", 8, "He$^+$", [("tb", 28)]),
    _r(("2body", 27), "H+", "O+", "chem", 8, "O", [("tb", 29)]),
    _r(("2body", 28), "O+", "H+", "chem", 8, "H", [("tb", 30)]),
    _r(("2body", 29), "O+", "OHx", "chem", 8, "H$_2$", [("tb", 3)]),
    _r(("2body", 30), "O+", "O", "chem", 8, "H$_2$", [("tb", 4)]),
    # CO-HCO+ cycle
    _r(("2body", 2), "CO", "HCO+", "chem", 9, "H$_3^+$", [("tb", 5)]),
    _r(("2body", 6), "C+", "HCO+", "chem", 9, "OH$_x$", [("tb", 10)]),
    _r(("2body", 6), "OHx", "HCO+", "chem", 9, "C$^+$", [("tb", 10)]),
    _r(("2body", 12), "HCO+", "CO", "rec", 9, despotic=[("tb", 17)]),
    # collisional
    _r((("2body", 15), ("2body", 16)), "H2", "H", "coll", 10, "H, H$_2$",
       [("tb", 21), ("tb", 22)]),
    _r(("2body", 17), "H", "H+", "coll", 10, "e", [("tb", 23)]),
]

N_ATHENAK = {"cr": 7, "photo": 6, "grain": 5, "2body": 31}
DESPOTIC_KEYS = ([("cr", i) for i in range(6)] + [("cr", "H"), ("cr", "H2")]
                 + [("ph", i) for i in range(6)] + [("tb", i) for i in range(36)])


def validate(reactions=REACTIONS):
    """Raise AssertionError unless the table covers AthenaK's network exactly."""
    gow = [r for r in reactions if r.network == "gow17"]
    ids = [i for r in gow for i in r.athenak]
    expected = {(g, n) for g, count in N_ATHENAK.items() for n in range(count)}
    assert set(ids) == expected, sorted(expected ^ set(ids))
    # 2-body (6) is drawn from both reactants; every other id appears once.
    dup = {i for i in ids if ids.count(i) > 1}
    assert dup == {("2body", 6)}, dup
    assert len(expected) == 49
    keys = [k for r in gow for k in r.despotic]
    assert set(keys) == set(DESPOTIC_KEYS), sorted(set(DESPOTIC_KEYS) ^ set(keys), key=str)
    for r in reactions:
        assert r.src in BY_NAME and r.dst in BY_NAME, r
        assert r.process in PROCESSES, r
        assert r.step in STEPS, r
    return True
