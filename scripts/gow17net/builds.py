"""The slide build sequence: one physical process per step, with a caption.

Steps and the reactions each introduces are in reactions.py (Reaction.step).
Captions here are also used by the movie.
"""
import pathlib

from . import render
from .reactions import STEPS

CAPTIONS = {
    1: "Twelve species are integrated; H, He, C, O, Si and the electrons follow "
       "from element and charge conservation.",
    2: "Cosmic rays reach deep into clouds, ionizing H, H$_2$, He, C and Si.",
    3: "Far-ultraviolet photons ionize C and Si and dissociate H$_2$, CO, CH$_x$ "
       "and OH$_x$; cosmic rays also make a weak ultraviolet field that breaks CO.",
    4: "Dust grains are where H$_2$ forms, and grain surfaces neutralize "
       "H$^+$, He$^+$, C$^+$ and Si$^+$.",
    5: "Free electrons recombine with the ions; H$_3^+$ recombines to H$_2$ + H "
       "or to three H atoms.",
    6: "H$_2^+$ reacts with H$_2$ to make H$_3^+$, which drives the ion chemistry; "
       "He$^+$ breaks H$_2$ apart.",
    7: "H$_3^+$ + C and C$^+$ + H$_2$ make the CH hydrides (CH$_x$), CH$_x$ + O makes "
       "CO, and He$^+$ + CO returns the carbon to C$^+$.",
    8: "Charge exchange with H$^+$ makes O$^+$; O$^+$ + H$_2$ and H$_3^+$ + O make the "
       "OH hydrides (OH$_x$), and OH$_x$ + C makes CO.",
    9: "CO and HCO$^+$ turn into each other: H$_3^+$ and cosmic rays make HCO$^+$ "
       "from CO, and HCO$^+$ recombines back to CO.",
    10: "In hot gas, collisions dissociate H$_2$ and ionize H.",
    11: "AthenaK updates the species in this order (block Gauss–Seidel); H$_2$ "
        "and CO limit the substep. Halos show how TIGRESS-NCR treats the species it keeps.",
    12: "TIGRESS-NCR adds EUV photoionization of H and H$_2$, and above 2$\\times$10$^4$ K "
        "takes the He and metal ionization state from CIE tables.",
}
assert set(CAPTIONS) == set(STEPS)


def render_sequence(layout, outdir, theme="light", widths=None):
    """Write build_NN_<name>.png for every step; return the paths."""
    outdir = pathlib.Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    paths = []
    for step, name in STEPS.items():
        fig, _ = render.draw(layout, step=step, widths=widths, theme=theme,
                             title=f"{step}. {name}", caption=CAPTIONS[step])
        path = outdir / f"build_{step:02d}_{name.replace(' ', '_').replace('+', 'p')}.png"
        render.save(fig, path)
        paths.append(path)
    return paths
