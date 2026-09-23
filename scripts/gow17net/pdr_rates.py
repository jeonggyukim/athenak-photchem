"""Per-reaction rates through the slab PDR of layouts/slab.py.

At each A_V of the slab profile the network is rebuilt with the slab's
parameters (n_H, T, chi, xi_CR, and N_H = 1.87e21 A_V for self-shielding) and
every reaction's rate per H nucleus [s^-1] is evaluated at the profile's
equilibrium abundances; the abundances are not re-solved. Rates are cached in
figures/gow17net/pdr_rates.npz.

A reaction destroys a species when its net stoichiometry for that species is
negative and forms it when positive, taken from DESPOTIC's reaction matrices
(`mat`). Net stoichiometry keeps a species that appears on both sides, such as
H2 in H2 + H2 -> H2 + 2H, from being counted as destroyed twice.
"""
import numpy as np

from . import fluxes
from .layouts import slab
from .pdr_profile import OUTDIR
from .reactions import REACTIONS

CACHE = OUTDIR / "pdr_rates.npz"
SPECIES = fluxes.gowmod.specListExtended

CLASSES = ("FUV", "CR", "grain", "recombination", "chemistry", "collisional")
_PROCESS_TO_CLASS = {"uv": "FUV", "cr": "CR", "crphoto": "CR", "grain": "grain",
                     "rec": "recombination", "chem": "chemistry",
                     "coll": "collisional"}


def _key_str(key):
    return f"{key[0]}:{key[1]}"


def _key_parse(s):
    g, i = s.split(":")
    return (g, int(i)) if i.isdigit() else (g, i)


def _compute(av, x):
    rates = {}
    for j, a in enumerate(av):
        net = fluxes.make_network(slab.NH, slab.T, slab.CHI, a, slab.XI_CR)
        xj = np.array([x[s][j] for s in fluxes.gowmod.specList])
        _, r = fluxes.reaction_rates(net, xj)
        for k, v in r.items():
            rates.setdefault(k, np.empty(len(av)))[j] = float(np.ravel(v)[0])
    return rates


def rates_vs_av():
    """(A_V array, {DESPOTIC key: rate per H nucleus [s^-1] at each A_V})."""
    av, x = slab.profile()
    params = [slab.NH, slab.T, slab.CHI, slab.XI_CR]
    if CACHE.exists():
        d = np.load(CACHE)
        if np.array_equal(d["AV"], av) and np.array_equal(d["params"], params):
            return av, {_key_parse(k): d["rates"][i] for i, k in enumerate(d["keys"])}
    rates = _compute(av, x)
    keys = list(rates)
    OUTDIR.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, AV=av, params=params, keys=[_key_str(k) for k in keys],
             rates=np.array([rates[k] for k in keys]))
    return av, rates


def net_stoich(key):
    """{species: net stoichiometry} of a DESPOTIC reaction (products > 0)."""
    group, i = key
    if group == "cr" and i == "H":
        return {"H": -1, "H+": 1, "e-": 1}
    if group == "cr" and i == "H2":
        return {"H2": -1, "H2+": 1, "e-": 1}
    gm = fluxes.gowmod
    rm = {"cr": gm._cr, "ph": gm._ph, "tb": gm._twobody}[group]
    col = np.asarray(rm.mat.todense() if hasattr(rm.mat, "todense") else rm.mat)[:, i]
    return {SPECIES[j]: int(col[j]) for j in range(len(SPECIES)) if col[j] != 0}


def _lhs(key):
    """{species: stoichiometry} on the left-hand side."""
    group, i = key
    if group == "cr" and i in ("H", "H2"):
        return {i: 1}
    gm = fluxes.gowmod
    rm = {"cr": gm._cr, "ph": gm._ph, "tb": gm._twobody}[group]
    return {SPECIES[j]: int(s) for j, s in zip(rm.lhsidx[i], rm.lhsstoich[i]) if s > 0}


def process_class(key):
    procs = {r.process for r in REACTIONS if r.network == "gow17" and key in r.despotic}
    assert len(procs) == 1, (key, procs)
    return _PROCESS_TO_CLASS[procs.pop()]


def _electron_stoich(key):
    """Electrons gained per reaction, from the change in total ion charge.

    DESPOTIC's lumped products CHx and OHx stand for ions that recombine, so
    reactions such as C+ + H2 -> CHx + H use an electron that the matrices
    omit; charge neutrality counts it.
    """
    return sum(v for s, v in net_stoich(key).items() if s.endswith("+"))


def by_reaction(species, sign):
    """{key: rate per H of `species` formed (sign=+1) or destroyed (sign=-1)}.

    For e- the stoichiometry is the change in total ion charge.
    """
    _, rates = rates_vs_av()
    out = {}
    for key, r in rates.items():
        st = (_electron_stoich(key) if species == "e-"
              else net_stoich(key).get(species, 0))
        if st * sign > 0:
            out[key] = abs(st) * r
    return out


def _by_process(species, sign):
    out = {}
    for key, r in by_reaction(species, sign).items():
        c = process_class(key)
        out[c] = out.get(c, 0.0) + r
    return out


def destruction_by_process(species):
    """{process class: destruction rate of `species` per H nucleus [s^-1]}."""
    return _by_process(species, -1)


def formation_by_process(species):
    """{process class: formation rate of `species` per H nucleus [s^-1]}."""
    return _by_process(species, +1)


def label(key):
    """Reaction as text, e.g. 'He+ + CO -> C+ + He + O'."""
    lhs = _lhs(key)
    rhs = dict(lhs)
    for s, v in net_stoich(key).items():
        rhs[s] = rhs.get(s, 0) + v
    lhs = [f"{v if v > 1 else ''}{s}" for s, v in lhs.items()]
    rhs = [f"{v if v > 1 else ''}{s}" for s, v in rhs.items() if v > 0]
    return " + ".join(lhs) + " -> " + " + ".join(rhs)


if __name__ == "__main__":
    av, rates = rates_vs_av()
    _, x = slab.profile()
    for s in ("H2", "C+", "C", "CO", "e-", "H+"):
        f = sum(formation_by_process(s).values())
        d = sum(destruction_by_process(s).values())
        print(f"{s:4s} max |F-D|/D = {np.max(np.abs(f - d) / d):.2e}")
