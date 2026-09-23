"""Colours, line styles and fonts shared by every layout."""

PROCESS_STYLE = {
    "cr": dict(color="#C0392B", ls="-", lw=1.7),
    "crphoto": dict(color="#C0392B", ls=(0, (4, 2)), lw=1.7),
    "uv": dict(color="#7D3C98", ls=(0, (4, 2)), lw=1.7),
    "grain": dict(color="#2E86C1", ls="-", lw=1.7),
    "rec": dict(color="#34495E", ls="-", lw=1.0),
    "chem": dict(color="#1B2631", ls="-", lw=1.5),
    "coll": dict(color="#E67E22", ls="-.", lw=1.5),
    "ncr": dict(color="#117A65", ls=(0, (6, 2, 1, 2)), lw=1.7),
}
PROCESS_LABEL = {
    "cr": "cosmic-ray ionization",
    "crphoto": "CR-induced photodissociation",
    "uv": "FUV photoreaction",
    "grain": "grain-assisted",
    "rec": "recombination (e)",
    "chem": "ion-molecule / neutral-neutral",
    "coll": "collisional (high T)",
    "ncr": "only in TIGRESS-NCR",
}

BOX_FILL = {"ion": "#D6EAF8", "molecule": "#D5F5E3", "ghost": "white"}
BOX_EDGE = {"ion": "#1B2631", "molecule": "#1B2631", "ghost": "#7F8C8D"}
NCR_HALO = {"non-equilibrium": "#F5B041", "ghost": "#AAB7B8", "equilibrium": "#48C9B0"}
NCR_LABEL = {"non-equilibrium": "NCR: non-equilibrium",
             "ghost": "NCR: ghost (element budget)",
             "equilibrium": "NCR: equilibrium or fixed ratio"}

# Build-sequence rendering: arrows introduced earlier than the current step.
PAST_COLOR = "#B3B6B7"
PARTNER_COLOR = "#5D6D7E"

FONT = dict(species=13, label=8, legend=11.5, caption=10.5, title=13)
