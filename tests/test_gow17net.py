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
