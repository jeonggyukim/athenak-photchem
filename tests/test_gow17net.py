import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "scripts"))

from gow17net import reactions, species  # noqa: E402


def test_reaction_table_covers_athenak_network():
    assert reactions.validate()


def test_every_drawn_species_is_used():
    used = {r.src for r in reactions.REACTIONS} | {r.dst for r in reactions.REACTIONS}
    assert set(species.DRAWN) <= used, set(species.DRAWN) - used
