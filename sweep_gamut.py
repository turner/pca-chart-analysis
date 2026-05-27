#!/usr/bin/env python3
"""Sweep (L*, ab_span) to find configurations that keep reference points in the sRGB gamut."""

import argparse
import numpy as np
import pandas as pd

_Xn, _Yn, _Zn = 0.95047, 1.00000, 1.08883
_M = np.array([
    [ 3.24062548, -1.53720797, -0.49862860],
    [-0.96893071,  1.87575606,  0.04151752],
    [ 0.05571012, -0.20402105,  1.05699594],
])

def _f_inv(t):
    delta = 6.0 / 29.0
    return np.where(t > delta, t ** 3, 3 * delta * delta * (t - 4.0 / 29.0))

def lab_to_linear_rgb_noclip(L, a, b):
    fy = (L + 16.0) / 116.0
    fx = a / 500.0 + fy
    fz = fy - b / 200.0
    X = _Xn * _f_inv(fx); Y = _Yn * _f_inv(fy); Z = _Zn * _f_inv(fz)
    XYZ = np.stack([X, Y, Z], axis=-1)
    return XYZ @ _M.T

def frac_in_gamut(x1, x2, L, ab_span, margin=0.1):
    x1m, x1M = x1.min() - margin, x1.max() + margin
    x2m, x2M = x2.min() - margin, x2.max() + margin
    a = ((x1 - x1m) / (x1M - x1m + 1e-12) - 0.5) * 2.0 * ab_span
    b = ((x2 - x2m) / (x2M - x2m + 1e-12) - 0.5) * 2.0 * ab_span
    L_arr = np.full_like(a, L)
    lin = lab_to_linear_rgb_noclip(L_arr, a, b)
    in_g = np.all((lin >= 0.0) & (lin <= 1.0), axis=-1)
    return float(in_g.mean())

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--target", type=float, default=0.95)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    x1 = df["x1"].to_numpy(float); x2 = df["x2"].to_numpy(float)

    Ls = list(range(40, 96, 5))
    spans = list(range(20, 121, 5))

    grid = np.zeros((len(Ls), len(spans)))
    for i, L in enumerate(Ls):
        for j, s in enumerate(spans):
            grid[i, j] = frac_in_gamut(x1, x2, L, s, args.margin)

    print("Fraction of 3122 reference haplotypes in sRGB gamut")
    print("rows = L*, cols = ab_span")
    print()
    header = "L* \\ span | " + " ".join(f"{s:4d}" for s in spans)
    print(header)
    print("-" * len(header))
    for i, L in enumerate(Ls):
        row = " ".join(f"{grid[i, j]*100:4.0f}" for j in range(len(spans)))
        print(f"   L*={L:3d}  | {row}")

    # find the configuration that meets target with the largest span;
    # among ties, prefer higher L* (pastel-ish, like the paper aims for).
    best = None
    for i, L in enumerate(Ls):
        for j, s in enumerate(spans):
            if grid[i, j] >= args.target:
                cand = (s, L)
                if best is None or cand > best[:2]:
                    best = (s, L, grid[i, j])
    if best is None:
        print(f"\nNo (L*, ab_span) pair reached target {args.target*100:.0f}% in this grid.")
    else:
        s, L, f = best
        print(f"\nBest @ ≥{args.target*100:.0f}% in-gamut: L*={L}, ab_span={s} -> {f*100:.1f}%")

    # also report: at each L*, the largest span that meets the target
    print(f"\nLargest ab_span meeting target ({args.target*100:.0f}%) at each L*:")
    for i, L in enumerate(Ls):
        ok = [spans[j] for j in range(len(spans)) if grid[i, j] >= args.target]
        if ok:
            print(f"  L*={L:3d}  max_span={max(ok):3d}  in-gamut={grid[i, spans.index(max(ok))]*100:.1f}%")
        else:
            print(f"  L*={L:3d}  (no span meets target)")

if __name__ == "__main__":
    main()
