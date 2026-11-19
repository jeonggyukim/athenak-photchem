import numpy as np, sys
sys.path.insert(0, ".")
from errors import read  # noqa
tr, ref = read("out_ref")
for run in ["out_f0.03_refresh_false", "out_f0.1_refresh_false"]:
    t, x = read(run)
    for k in ["CHx", "H3+", "OHx", "O+"]:
        e = x[k] / np.interp(t, tr, ref[k]) - 1
        i = np.nanargmax(np.abs(e[1:])) + 1
        print(run, k, f"t={t[i]:.4g} err={e[i]:.3g} x={x[k][i]:.3g} ref={np.interp(t[i], tr, ref[k]):.3g}",
              f"cycle={i}")
