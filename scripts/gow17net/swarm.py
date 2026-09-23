"""Atom swarm: carbon atoms hopping between species at each depth of the slab.

Each carbon atom is in one of C+, C, CHx, CO, HCO+. A reaction that moves
carbon from species s to t (net stoichiometry -1 on s, +1 on t) with rate R
per H nucleus moves each carbon atom in s at rate R / x_s. The resulting
continuous-time Markov chain has the equilibrium abundances as its
stationary distribution, so a swarm started all as C+ relaxes into the PDR
layers: fast near the surface, slowly at depth.
"""
import numpy as np
from scipy.linalg import expm

from . import pdr_rates
from .layouts import slab

CARBON = ["C+", "C", "CHx", "CO", "HCO+"]
YR = 3.156e7


def generator():
    """(A_V, Q) with Q[k] the 5x5 rate matrix [s^-1] at depth k (rows sum to 0)."""
    av, x = slab.profile()
    _, rates = pdr_rates.rates_vs_av()
    n = len(CARBON)
    Q = np.zeros((len(av), n, n))
    for key, r in rates.items():
        st = pdr_rates.net_stoich(key)
        src = [s for s in CARBON if st.get(s, 0) < 0]
        dst = [s for s in CARBON if st.get(s, 0) > 0]
        if len(src) != 1 or len(dst) != 1:
            continue      # reaction does not move carbon between two species
        i, j = CARBON.index(src[0]), CARBON.index(dst[0])
        xs = np.maximum(x[src[0]], 1e-40)
        Q[:, i, j] += r / xs
    for k in range(len(av)):
        np.fill_diagonal(Q[k], 0.0)
        Q[k][np.diag_indices(n)] = -Q[k].sum(axis=1)
    return av, Q


def stationary(Qk):
    """Stationary distribution of one generator matrix."""
    n = Qk.shape[0]
    A = np.vstack([Qk.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    return np.linalg.lstsq(A, b, rcond=None)[0]


def simulate(n_atoms=60, times_yr=None, seed=1, start="C+"):
    """States [frame, depth, atom] as indices into CARBON, at `times_yr`."""
    av, Q = generator()
    if times_yr is None:
        times_yr = np.logspace(0, 7, 140)
    rng = np.random.default_rng(seed)
    state = np.full((len(av), n_atoms), CARBON.index(start))
    out = np.empty((len(times_yr), len(av), n_atoms), dtype=np.int8)
    t_prev = 0.0
    for f, t in enumerate(times_yr):
        dt = (t - t_prev) * YR
        for k in range(len(av)):
            P = expm(Q[k] * dt)
            P = np.clip(P, 0.0, None)
            P /= P.sum(axis=1, keepdims=True)
            cum = np.cumsum(P[state[k]], axis=1)
            u = rng.random(n_atoms)[:, None]
            state[k] = (u > cum).sum(axis=1)
        out[f] = state
        t_prev = t
    return av, np.asarray(times_yr), out
