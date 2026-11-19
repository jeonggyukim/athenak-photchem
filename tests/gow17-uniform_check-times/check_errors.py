"""Errors at the check times of run.sh: for each condition, time and solver, the
floored error |x - x_ref| / max(x_ref, 1e-3 x_elem) per species (x_elem: element
abundance of the species' heaviest element, H 1, He 0.1, C 1.6e-4, O 3.2e-4,
Si 1.7e-6) and |e/e_ref - 1| for the internal energy. Run in OUT."""
import glob
import os
import numpy as np

ELEM = {"He+": 0.1, "OHx": 3.2e-4, "CHx": 1.6e-4, "CO": 1.6e-4, "C+": 1.6e-4,
        "HCO+": 1.6e-4, "H2": 1.0, "H+": 1.0, "H3+": 1.0, "H2+": 1.0, "O+": 3.2e-4,
        "Si+": 1.7e-6}


def final(d):
    f = sorted(glob.glob(os.path.join(d, "tab", "*.tab")))[-1]
    with open(f) as fh:
        t = float(fh.readline().split("time=")[1].split()[0])
        cols = fh.readline().split()[1:]
        row = [float(v) for v in fh.readline().split()]
    return t, {c.split("chem_")[-1]: v for c, v in zip(cols, row)}


def errors(run, ref):
    e = {k: abs(run[k] - ref[k]) / max(ref[k], 1e-3 * x) for k, x in ELEM.items()}
    e["eint"] = abs(run["eint"] / ref["eint"] - 1)
    return e


def main():
    tags = sorted({d.rsplit("_", 1)[0] for d in os.listdir("out")},
                  key=lambda s: (s.split("_t")[0], float(s.split("_t")[1])))
    print(f"{'case':14s} {'solver':6s} {'worst':>7s} {'species':>7s} {'eint':>7s}")
    for tag in tags:
        t, ref = final(f"out/{tag}_ref")
        for s in ["bdf", "si"]:
            t2, x = final(f"out/{tag}_{s}")
            assert abs(t2 - t) < 1e-12 * max(t, 1)
            e = errors(x, ref)
            k = max((k for k in e if k != "eint"), key=e.get)
            print(f"{tag:14s} {s:6s} {e[k]:7.3f} {k:>7s} {e['eint']:7.3f}")


if __name__ == "__main__":
    main()
