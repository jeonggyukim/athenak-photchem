import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

from gow17net import arcs, reactions, species  # noqa: E402


def test_reaction_table_covers_athenak_network():
    assert reactions.validate()


def test_every_drawn_species_is_used():
    used = {r.src for r in reactions.REACTIONS} | {r.dst for r in reactions.REACTIONS}
    assert set(species.DRAWN) <= used, set(species.DRAWN) - used


def test_parallel_arrows_bend_apart():
    rads = arcs.assign(reactions.REACTIONS)
    groups = {}
    for r, rad in zip(reactions.REACTIONS, rads):
        pair, forward = arcs.canonical(r.src, r.dst)
        # Express every rad in the canonical direction before comparing.
        groups.setdefault(pair, []).append(rad if forward else -rad)
    for pair, offs in groups.items():
        assert len(set(offs)) == len(offs), (pair, offs)


def test_offsets_symmetric():
    for n in range(1, 6):
        o = arcs.offsets(n, 0.25)
        assert len(set(o)) == n
        assert abs(sum(o)) < 1e-12


LAYOUTS = ("grid", "grid_charge_rows", "radial", "slab", "ladder")


def _layout(name):
    from gow17net import layouts
    return layouts.get(name)


def test_boxes_do_not_overlap_and_fit():
    for name in LAYOUTS:
        lay = _layout(name)
        w, h = lay.box_w, lay.box_h
        pos = lay.positions
        assert set(pos) == set(species.DRAWN), (name, set(species.DRAWN) ^ set(pos))
        names = list(pos)
        for i, a in enumerate(names):
            xa, ya = pos[a]
            assert lay.xlim[0] + w / 2 <= xa <= lay.xlim[1] - w / 2, (name, a)
            assert lay.ylim[0] + h / 2 <= ya <= lay.ylim[1] - h / 2, (name, a)
            for b in names[i + 1:]:
                xb, yb = pos[b]
                apart = abs(xa - xb) >= w + 0.1 or abs(ya - yb) >= h + 0.1
                assert apart, (name, a, b)


def test_boxes_stay_put_across_build_steps():
    from gow17net import render
    lay = _layout("grid")
    extents = None
    for step in (None, 1, 5, 12):
        fig, boxes = render.draw(lay, step=step)
        fig.canvas.draw()
        ext = {k: tuple(round(v, 3) for v in b.get_window_extent().bounds)
               for k, b in boxes.items()}
        import matplotlib.pyplot as plt
        plt.close(fig)
        if extents is None:
            extents = ext
        assert ext == extents, step


def test_width_fractions_are_shares():
    from gow17net import fluxes, widths
    for phase in fluxes.PHASES:
        frac = widths.phase_fractions(phase)
        assert frac, phase
        for i, f in frac.items():
            assert 0.0 <= f <= 1.0 + 1e-9, (phase, i, f)
