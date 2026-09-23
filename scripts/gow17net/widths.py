"""Arrow widths from reaction rates at equilibrium.

An arrow's width is its reaction's largest share of destruction among its
reactants (electrons excluded). Using only the drawn source hides dominant
sinks drawn from the other reactant: in molecular gas He+ + CO -> C+ carries
most of He+'s destruction but is drawn from CO, and H3+ + CO -> HCO+ is drawn
from CO too. With the largest share, such arrows come out thick.
"""
from . import fluxes
from .reactions import REACTIONS

LW_MIN, LW_MAX = 0.35, 6.0


def _reactants(key):
    """{species: stoichiometry} on the left-hand side of a DESPOTIC reaction."""
    group, i = key
    if group == "cr" and i in ("H", "H2"):
        return {i: 1}
    gm = fluxes.gowmod
    rm = {"cr": gm._cr, "ph": gm._ph, "tb": gm._twobody}[group]
    names = gm.specListExtended
    return {names[j]: int(s) for j, s in zip(rm.lhsidx[i], rm.lhsstoich[i]) if s > 0}


def fractions(rates):
    """{reaction index: branching fraction} from DESPOTIC per-reaction rates."""
    destruction = {}
    for key, rate in rates.items():
        for sp, st in _reactants(key).items():
            destruction[sp] = destruction.get(sp, 0.0) + st * rate
    out = {}
    for i, r in enumerate(REACTIONS):
        if r.network != "gow17":
            continue
        shares = {}
        for k in r.despotic:
            for sp, st in _reactants(k).items():
                if sp != "e-" and destruction.get(sp, 0.0) > 0.0:
                    shares[sp] = shares.get(sp, 0.0) + st * rates[k] / destruction[sp]
        out[i] = max(shares.values(), default=0.0)
    return out


def phase_fractions(phase):
    _, rates = fluxes.phase_rates(phase)
    return fractions(rates)


def to_widths(frac, lw_min=LW_MIN, lw_max=LW_MAX):
    return {i: lw_min + (lw_max - lw_min) * min(f, 1.0) for i, f in frac.items()}
