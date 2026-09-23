"""Equilibrium abundances of the GOW17 network through a 1D PDR.

A slab with fixed n_H = 1e3 cm^-3 and T = 50 K, lit by an incident FUV field
chi = 1 and a cosmic-ray ionization rate xi_CR = 2e-16 s^-1 per H. Temperature
is held fixed rather than solved from thermal balance, so the profile isolates
the chemistry. Each A_V point is an independent equilibrium: the FUV rates are
attenuated by exp(-gamma A_V) and the self-shielding columns use
N_H = 1.87e21 A_V cm^-2.

The profile is cached in figures/gow17net/pdr_profile.npz. Run as a script to
write figures/gow17net/pdr_profile.png.
"""
import pathlib

import numpy as np

from . import fluxes

NH = 1.0e3
T = 50.0
CHI = 1.0
XI_CR = 2.0e-16
AV = np.logspace(-2, 1, 40)

OUTDIR = pathlib.Path(__file__).resolve().parents[2] / "figures" / "gow17net"
CACHE = OUTDIR / "pdr_profile.npz"
SPECIES = fluxes.gowmod.specListExtended


def _compute():
    x = np.empty((len(SPECIES), len(AV)))
    resid = np.empty(len(AV))
    for j, av in enumerate(AV):
        net = fluxes.make_network(NH, T, CHI, av, XI_CR)
        xj = fluxes.equilibrium(net)
        x[:, j] = net.extendAbundances(xj, outdict=True).x.ravel()
        # Largest |dx/dt| * t_eq / x over species above 1e-12, with t_eq = 1 Myr.
        dxdt = np.asarray(net.dxdt(xj, 0.0)).ravel()
        ok = xj > 1.0e-12
        resid[j] = np.max(np.abs(dxdt[ok]) * 3.15e13 / xj[ok])
    return x, resid


def _load():
    if CACHE.exists():
        d = np.load(CACHE)
        if np.array_equal(d["AV"], AV) and list(d["species"]) == SPECIES \
                and np.array_equal(d["params"], [NH, T, CHI, XI_CR]):
            return d["x"], d["resid"]
    x, resid = _compute()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    np.savez(CACHE, AV=AV, species=SPECIES, x=x, resid=resid,
             params=[NH, T, CHI, XI_CR])
    return x, resid


def profile():
    """(A_V array, dict species -> equilibrium abundance per H nucleus)."""
    x, _ = _load()
    return AV.copy(), {s: x[i] for i, s in enumerate(SPECIES)}


def residual():
    """Per A_V, the largest fractional change of any species over 1 Myr at the
    returned state; values << 1 mean the equilibrium has converged."""
    return _load()[1]


def peak_av():
    """dict species -> A_V at which its abundance is largest."""
    av, x = profile()
    return {s: av[np.argmax(v)] for s, v in x.items()}


def _crossing(av, a, b):
    """First A_V where a - b changes sign, interpolated in log A_V and log(a/b);
    NaN if the sign never changes on the grid."""
    r = np.log(np.maximum(a, 1e-300)) - np.log(np.maximum(b, 1e-300))
    idx = np.nonzero(np.sign(r[:-1]) != np.sign(r[1:]))[0]
    if len(idx) == 0:
        return np.nan
    i = idx[0]
    la = np.log10(av[i:i + 2])
    return 10 ** (la[0] + (la[1] - la[0]) * r[i] / (r[i] - r[i + 1]))


def transitions():
    """A_V of the H/H2, C+/C and C/CO (or C+/CO) transitions."""
    av, x = profile()
    out = {"H/H2": _crossing(av, x["H"], 2 * x["H2"]),
           "C+/C": _crossing(av, x["C+"], x["C"])}
    c_dominates = np.any((x["C"] > x["C+"]) & (x["C"] > x["CO"]))
    if c_dominates:
        out["C/CO"] = _crossing(av, x["C"], x["CO"])
    else:
        out["C+/CO"] = _crossing(av, x["C+"], x["CO"])
    return out


PLOT = [("H", "H", 1), ("H2", "2 H$_2$", 2), ("C+", "C$^+$", 1), ("C", "C", 1),
        ("CO", "CO", 1), ("e-", "e$^-$", 1), ("OHx", "OH$_x$", 1),
        ("HCO+", "HCO$^+$", 1)]


def plot(path=OUTDIR / "pdr_profile.png"):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    av, x = profile()
    fig, ax = plt.subplots(figsize=(7, 5))
    for s, label, fac in PLOT:
        ax.loglog(av, fac * x[s], label=label, lw=1.8)
    for k, (name, a) in enumerate(transitions().items()):
        if np.isfinite(a):
            ax.axvline(a, color="0.4", ls="--", lw=1)
            ax.text(a * 1.05, 10 ** (-1.5 - 1.2 * k), f"{name}\n$A_V$={a:.2f}",
                    fontsize=9, color="0.25")
        else:
            ax.text(av[0] * 1.1, 10 ** (-1.5 - 1.2 * k),
                    f"{name}: $A_V$ < {av[0]:g}", fontsize=9, color="0.25")
    ax.set_xlim(av[0], av[-1])
    ax.set_ylim(1e-12, 2)
    ax.set_xlabel("$A_V$ [mag]")
    ax.set_ylabel("abundance per H nucleus")
    ax.set_title(rf"GOW17 equilibrium: $n_{{\rm H}}$={NH:.0e} cm$^{{-3}}$, "
                 rf"$T$={T:.0f} K, $\chi$={CHI:g}, $\xi_{{\rm CR}}$={XI_CR:.0e} s$^{{-1}}$",
                 fontsize=10)
    ax.legend(ncol=2, fontsize=9, loc="lower left")
    fig.tight_layout()
    OUTDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


if __name__ == "__main__":
    print(plot())
    print("transitions:", {k: round(float(v), 3) for k, v in transitions().items()})
    pk = peak_av()
    print("peaks:", ", ".join(f"{s}={a:.3g}" for s, a in sorted(pk.items(), key=lambda kv: kv[1])))
    print("max residual:", residual().max())
