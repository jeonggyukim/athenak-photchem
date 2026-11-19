"""Worst relative error over 30 Myr per species and eint, refresh_rates true vs false,
(optional argument: run prefix, e.g. nocap_)
against out_ref. Usage: run in ~/Documents/Athenak-chem-tests/gow17-refresh-rates/ (the run output)."""
import glob
import numpy as np


def read(d):
    t, rows, cols = [], [], None
    for f in sorted(glob.glob(f"{d}/tab/*.tab")):
        with open(f) as fh:
            head = fh.readline()
            cols = fh.readline().split()[1:]
            rows.append([float(v) for v in fh.readline().split()])
        t.append(float(head.split("time=")[1].split()[0]))
    a = np.array(rows)
    names = [c.split("chem_")[-1] for c in cols]
    return np.array(t), {n: a[:, i] for i, n in enumerate(names)}


import sys
PREFIX = sys.argv[1] if len(sys.argv) > 1 else ""
tr, ref = read("out_ref")
keys = ["eint"] + [k for k in ref if k not in ("gid", "i", "x1v", "dens", "velx", "vely",
                                               "velz", "eint")]
print("run".ljust(22), " ".join(k.rjust(7) for k in keys), "  worst  final_worst")
for f in ["0.3", "0.1", "0.03"]:
    for r in ["false", "true"]:
        t, x = read(f"out_{PREFIX}f{f}_refresh_{r}")
        m = t > 0
        e = {k: np.abs(x[k][m] / np.interp(t[m], tr, ref[k]) - 1) for k in keys}
        fin = {k: abs(x[k][-1] / ref[k][-1] - 1) for k in keys}
        print(f"f{f} refresh {r}".ljust(22),
              " ".join(f"{np.max(e[k]):7.3f}" for k in keys),
              f"  {max(np.max(v) for v in e.values()):.3f} {max(fin.values()):.4f}")
