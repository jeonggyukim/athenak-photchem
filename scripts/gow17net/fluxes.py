"""Per-reaction fluxes of the GOW17 network at chemical equilibrium.

Uses the pure-numpy GOW17 implementation in DESPOTIC (Krumholz), which carries
the same 12 integrated and 6 ghost species as AthenaK. For each phase the
species are integrated to steady state at fixed temperature, and every
reaction's rate per H nucleus [s^-1] is returned, in DESPOTIC's own numbering.

DESPOTIC must be importable; set DESPOTIC_ROOT if it is not at
~/Dropbox/Projects/despotic.
"""
import contextlib
import os
import pathlib
import sys
import warnings

import numpy as np
from scipy.integrate import solve_ivp

DESPOTIC_ROOT = pathlib.Path(
    os.environ.get("DESPOTIC_ROOT", pathlib.Path.home() / "Dropbox/Projects/despotic"))
sys.path.insert(0, str(DESPOTIC_ROOT))


@contextlib.contextmanager
def _quiet_fds():
    """Send fd 1 and 2 to /dev/null, including output of child processes."""
    sys.stdout.flush()
    sys.stderr.flush()
    saved = [os.dup(1), os.dup(2)]
    devnull = os.open(os.devnull, os.O_WRONLY)
    try:
        os.dup2(devnull, 1)
        os.dup2(devnull, 2)
        yield
    finally:
        os.dup2(saved[0], 1)
        os.dup2(saved[1], 2)
        for fd in saved + [devnull]:
            os.close(fd)


# Importing despotic tries to build an unrelated C++ extension (pwind), which
# fails without GSL and prints the compiler errors; the GOW network is unaffected.
with warnings.catch_warnings(), _quiet_fds():
    warnings.simplefilter("ignore")
    import despotic.chemistry.GOW  # noqa: E402,F401
# despotic.chemistry re-exports the GOW class under the module's own name, so
# the attribute path yields the class; the module itself is in sys.modules.
gowmod = sys.modules["despotic.chemistry.GOW"]


def _scalar(f):
    return lambda *a, **k: float(np.ravel(f(*a, **k))[0])


# DESPOTIC's shielding functions return 1-element arrays, which numpy 2 will
# not store into a float array slot; wrap them to return scalars.
gowmod._ph.shield_func = [_scalar(f) for f in gowmod._ph.shield_func]

KB = gowmod.kB
MH = gowmod.mH
AV_TO_NH = 1.87e21  # cm^-2 per magnitude, as GOW17 uses

# Canonical ISM phases. The FUV field is attenuated by exp(-av_fac * AV) in the
# photorates and by the column NH in the self-shielding functions.
PHASES = {
    "WNM": dict(nH=0.5, T=8000.0, chi=1.0, AV=0.0, ionrate=2.0e-16),
    "CNM": dict(nH=30.0, T=80.0, chi=1.0, AV=0.0, ionrate=2.0e-16),
    "molecular": dict(nH=1.0e3, T=15.0, chi=1.0, AV=5.0, ionrate=1.0e-17),
}
IONRATE = 2.0e-16


def make_network(nH, T, chi, AV, ionrate=IONRATE):
    net = gowmod.GOW(info={"Tg": T, "AV": AV})
    net.nH = nH
    net.temp = T
    net.chi = chi
    net.ionRate = ionrate
    net.NH = max(AV * AV_TO_NH, gowmod._small)
    return net


def equilibrium(net, t_end=3.0e15):
    """Integrate the 12 species to steady state at fixed temperature."""
    spec = gowmod.specList
    x0 = np.full(len(spec), 1.0e-12)
    x0[spec.index("C+")] = net.xC
    x0[spec.index("Si+")] = net.xSi
    x0[spec.index("H+")] = 1.0e-4
    sol = solve_ivp(lambda t, x: net.dxdt(x, t), (0.0, t_end), x0,
                    method="BDF", rtol=1.0e-6, atol=1.0e-25)
    return np.clip(sol.y[:, -1], 0.0, None)


def reaction_rates(net, x):
    """Rate per H nucleus [s^-1] of every reaction, keyed by (group, index).

    Mirrors GOW.dxdt, but keeps each reaction's rate instead of summing them
    through the stoichiometric matrix. The H and H2 cosmic-ray ionizations are
    handled outside DESPOTIC's reaction tables and are returned as
    ("cr", "H") and ("cr", "H2").
    """
    abd = net.extendAbundances(x, outdict=True)
    xe = abd.x
    n = net.nH * net.cfac
    rates = {}

    def lhs_product(rm, coef):
        return coef * np.prod(xe[rm.lhsidx] ** rm.lhsstoich, axis=1).ravel() \
            * n ** (rm.lhssum - 1)

    cr = gowmod._cr
    for i, r in enumerate(lhs_product(cr, cr.rate_fac * net.ionRate)):
        rates[("cr", i)] = r
    gamma_h = (2.3 * abd["H2"] + 1.5 * abd["H"]) * net.ionRate
    rates[("cr", "H")] = gamma_h * abd["H"]
    rates[("cr", "H2")] = 2.0 * gamma_h * abd["H2"]

    ph = gowmod._ph
    sigma_tot = np.sqrt(net.sigmaNT ** 2 + KB * net.temp / (net.mu(x) * MH))
    shield_args = [[net.NH * abd["C"], net.NH * abd["H2"]],
                   [net.NH * abd["CO"], net.NH * abd["H2"]],
                   [net.NH * abd["H2"], sigma_tot]]
    coef = ph.rate_fac * net.chi * np.exp(-net.AV * ph.av_fac)
    if len(ph.shield_func) > 0:
        fshield = np.array([f(*a) for f, a in zip(ph.shield_func, shield_args)])
        coef[ph.shield_func_idx] *= fshield
    for i, r in enumerate(lhs_product(ph, coef)):
        rates[("ph", i)] = r

    tb = gowmod._twobody
    k = gowmod._twobody_ratecoef(net.temp, n, abd["H"] * n, abd["H2"] * n,
                                 abd["e-"] * n,
                                 np.exp(-1.87 * net.AV) * net.chi)
    kgas = np.where(tb.grain, 0.0, k)
    gas = lhs_product(tb, kgas)
    grain = np.where(tb.grain, k * net.Zd, 0.0) * xe[tb.lhsidx[:, 0]].ravel() * n
    for i, r in enumerate(gas + grain):
        rates[("tb", i)] = r
    return abd, rates


def phase_rates(name):
    p = PHASES[name]
    net = make_network(p["nH"], p["T"], p["chi"], p["AV"],
                       p.get("ionrate", IONRATE))
    x = equilibrium(net)
    return reaction_rates(net, x)
