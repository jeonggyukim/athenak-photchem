#!/usr/bin/env python
"""Render GOW17 network diagrams with the gow17net package.

Examples:
    gow17net_cli.py --layout grid                      # full network
    gow17net_cli.py --layout grid --builds             # 12-slide build sequence
    gow17net_cli.py --layout grid --phase WNM CNM molecular
    gow17net_cli.py --layout radial --theme dark --step 7
"""
import argparse
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from gow17net import builds, layouts, render, widths  # noqa: E402

FIGDIR = pathlib.Path(__file__).resolve().parents[1] / "figures" / "gow17net"


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layout", default="grid")
    p.add_argument("--phase", nargs="+", help="arrow widths from these phases")
    p.add_argument("--step", type=int, help="one build step")
    p.add_argument("--builds", action="store_true", help="every build step")
    p.add_argument("--theme", default="light", choices=list(render.THEMES))
    p.add_argument("--outdir", type=pathlib.Path)
    args = p.parse_args()

    layout = layouts.get(args.layout)
    outdir = args.outdir or FIGDIR / args.layout
    outdir.mkdir(parents=True, exist_ok=True)
    suffix = "" if args.theme == "light" else f"_{args.theme}"

    if args.builds:
        for path in builds.render_sequence(layout, outdir / f"builds{suffix}",
                                           theme=args.theme):
            print("wrote", path)
        return
    if args.phase:
        for phase in args.phase:
            w = widths.to_widths(widths.phase_fractions(phase))
            fig, _ = render.draw(layout, widths=w, theme=args.theme, step=args.step,
                                 title=f"GOW17 network, {phase}: arrow width = largest "
                                       "share of a reactant's destruction")
            path = outdir / f"network_{phase}{suffix}.png"
            render.save(fig, path)
            print("wrote", path)
        return
    fig, _ = render.draw(layout, step=args.step, theme=args.theme,
                         caption=builds.CAPTIONS.get(args.step))
    name = f"network_step{args.step:02d}" if args.step else "network"
    path = outdir / f"{name}{suffix}.png"
    render.save(fig, path)
    print("wrote", path)


if __name__ == "__main__":
    main()
