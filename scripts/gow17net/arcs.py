"""Arc offsets for arrows that join the same two species.

matplotlib's arc3 bends an arrow to the right of its direction of travel for
rad > 0. Two arrows between A and B drawn in opposite directions with rads of
opposite sign therefore bend to the same side and overlap. Offsets are assigned
in a canonical direction (from the species that comes first in SPECIES) and
negated for arrows drawn the other way.
"""
from .reactions import PROCESSES
from .species import SPECIES

_ORDER = {s.name: i for i, s in enumerate(SPECIES)}
_PROCESS_ORDER = {p: i for i, p in enumerate(PROCESSES)}


def canonical(src, dst):
    """(first, second) of the pair in SPECIES order, and whether src is first."""
    if _ORDER[src] <= _ORDER[dst]:
        return (src, dst), True
    return (dst, src), False


def offsets(n, step):
    """n distinct offsets, symmetric about zero; zero included when n is odd."""
    if n == 1:
        return [0.0]
    out = [0.0] if n % 2 else []
    k = 1
    while len(out) < n:
        mag = step * (k if n % 2 else k - 0.5)
        out += [mag, -mag]
        k += 1
    return out[:n]


def assign(reactions, step=0.25, overrides=None):
    """rad for every reaction, as a list parallel to `reactions`.

    `overrides` maps an index in `reactions` to a rad (in the arrow's own
    direction) and is applied last.
    """
    groups = {}
    for i, r in enumerate(reactions):
        pair, _ = canonical(r.src, r.dst)
        groups.setdefault(pair, []).append(i)
    rads = [0.0] * len(reactions)
    for idx in groups.values():
        idx.sort(key=lambda i: (_PROCESS_ORDER.get(reactions[i].process, 99),
                                reactions[i].athenak))
        for i, off in zip(idx, offsets(len(idx), step)):
            _, forward = canonical(reactions[i].src, reactions[i].dst)
            rads[i] = off if forward else -off
    for i, rad in (overrides or {}).items():
        rads[i] = rad
    return rads
