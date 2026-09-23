"""Write figures/gow17net/index.html: every figure and movie, with a line each.

    python -m gow17net.contact_sheet     (from scripts/)
"""
import html
import pathlib

FIG = pathlib.Path(__file__).resolve().parents[2] / "figures" / "gow17net"

ENTRIES = [
    ("Reaction network (grid layout)", [
        ("grid/network.png", "Full GOW17 network, 49 reactions, solver annotations."),
        ("grid/network_step07_dark.png", "Dark theme, build step 7 (carbon chemistry)."),
        ("grid/network_CNM.png", "Arrow width = largest share of a reactant's destruction, CNM."),
        ("check/grid_molecular.png", "Same, molecular phase."),
        ("grid/build_01_species.png", "Build sequence, step 1 of 12 (see grid/ for all)."),
        ("grid/build_07_carbon_chemistry.png", "Build step 7."),
    ]),
    ("Other layouts", [
        ("check/radial_full.png", "Spherical PDR prototype (crowded; not chosen)."),
        ("check/slab/slab_full.png", "Plane-parallel slab with boxes (dropped after feedback)."),
        ("check/ladder/ladder_step12.png", "Ionization ladder for warm/hot gas, NCR-only step."),
    ]),
    ("PDR structure", [
        ("check/streams/streams_v2.png", "Element streams: share of each element per species vs A_V."),
        ("check/pdr_why/fig_dominant_process.png",
         "Which process controls each species at each depth (FUV to chemistry near A_V ~ 1)."),
        ("check/pdr_why/fig_rate_crossings.png", "Reaction rates behind it, for H2, C+ and CO."),
        ("pdr_profile.png", "One-zone equilibrium profile, n_H = 1e3, chi = 1."),
        ("check/pdr1d/compare_onezone.png", "1D march vs one-zone (task R1)."),
    ]),
    ("Movies (1080p)", [
        ("movie/BuildNetwork.mp4", "Network built one process at a time, grid layout."),
        ("movie/slab/BuildNetwork.mp4", "Same on the plane-parallel slab (superseded)."),
        ("movie/swarm.mp4", "Carbon atoms settling into the PDR layers, 1 yr to 2 Myr."),
    ]),
]


def main():
    parts = ["<!doctype html><meta charset=utf-8><title>gow17net figures</title>",
             "<style>body{font-family:system-ui;margin:24px;background:#fafafa}"
             "h2{margin-top:32px}.grid{display:flex;flex-wrap:wrap;gap:18px}"
             ".card{width:420px;background:white;padding:10px;border:1px solid #ddd}"
             ".card img,.card video{width:100%}.card p{font-size:14px;margin:6px 0 0}"
             ".missing{color:#999}</style>",
             "<h1>GOW17 network figures and movies</h1>"]
    for section, items in ENTRIES:
        parts.append(f"<h2>{html.escape(section)}</h2><div class=grid>")
        for rel, text in items:
            path = FIG / rel
            if not path.exists():
                parts.append(f"<div class='card missing'><p>{html.escape(rel)}: not made yet. "
                             f"{html.escape(text)}</p></div>")
                continue
            media = (f"<video src='{rel}' controls muted></video>" if rel.endswith(".mp4")
                     else f"<a href='{rel}'><img src='{rel}' loading=lazy></a>")
            parts.append(f"<div class=card>{media}<p><b>{html.escape(rel)}</b><br>"
                         f"{html.escape(text)}</p></div>")
        parts.append("</div>")
    out = FIG / "index.html"
    out.write_text("\n".join(parts))
    print("wrote", out)


if __name__ == "__main__":
    main()
