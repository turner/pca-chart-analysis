#!/usr/bin/env python3
"""
Binary gamut-membership map in PCA (x1, x2) coordinates.

For each pixel of the PCA bounding box:
  1. Apply the paper's formula:
       a* = ((x1 - x1_min)/(x1_max - x1_min) - 0.5) * 2 * ab_span
       b* = ((x2 - x2_min)/(x2_max - x2_min) - 0.5) * 2 * ab_span
       L* = L_fixed
  2. Convert Lab -> linear sRGB WITHOUT clipping (manual D65 pipeline).
  3. If all three linear-RGB channels are in [0, 1]: paint WHITE (in gamut).
     Otherwise: paint PURE RED (was out of gamut, would be clipped).

Input: a TSV with x1, x2 columns (e.g. reference_pca_metadata.tsv).
Output: <out>.png and <out>.pdf.
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# Lab -> linear sRGB, no clipping (D65, IEC 61966-2-1) ---------------------
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
    X = _Xn * _f_inv(fx)
    Y = _Yn * _f_inv(fy)
    Z = _Zn * _f_inv(fz)
    XYZ = np.stack([X, Y, Z], axis=-1)
    return XYZ @ _M.T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", default="gamut_mask")
    ap.add_argument("--L", type=float, default=80.0)
    ap.add_argument("--ab-span", type=float, default=120.0)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--grid", type=int, default=800)
    ap.add_argument("--figwidth", type=float, default=7.5)
    ap.add_argument("--figheight", type=float, default=6.5)
    ap.add_argument("--dpi", type=int, default=180)
    ap.add_argument("--show-points", action="store_true",
                    help="Overlay reference haplotypes as small black dots")
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    if "x1" not in df.columns or "x2" not in df.columns:
        raise SystemExit(f"TSV must contain columns x1 and x2; got {list(df.columns)}")
    x1 = df["x1"].to_numpy(float)
    x2 = df["x2"].to_numpy(float)

    x1_min, x1_max = x1.min() - args.margin, x1.max() + args.margin
    x2_min, x2_max = x2.min() - args.margin, x2.max() + args.margin

    gx = np.linspace(x1_min, x1_max, args.grid)
    gy = np.linspace(x2_min, x2_max, args.grid)
    GX, GY = np.meshgrid(gx, gy)

    A = ((GX - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    B = ((GY - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    L = np.full_like(A, args.L)

    lin = lab_to_linear_rgb_noclip(L, A, B)
    in_gamut = np.all((lin >= 0.0) & (lin <= 1.0), axis=-1)

    # white where in gamut, pure red where out of gamut
    img = np.zeros((args.grid, args.grid, 3), dtype=float)
    img[in_gamut] = [1.0, 1.0, 1.0]
    img[~in_gamut] = [1.0, 0.0, 0.0]

    fig, ax = plt.subplots(figsize=(args.figwidth, args.figheight), dpi=args.dpi)
    ax.imshow(
        img,
        extent=(x1_min, x1_max, x2_min, x2_max),
        origin="lower",
        aspect="auto",
        interpolation="nearest",
    )

    if args.show_points:
        ax.scatter(x1, x2, c="black", s=2.0, alpha=0.6, linewidths=0)

    frac_in = float(in_gamut.mean())
    ax.set_xlim(x1_min, x1_max)
    ax.set_ylim(x2_min, x2_max)

    ax.set_xlabel("PC1")
    ax.set_ylabel("PC2")
    ax.set_title(
        f"sRGB gamut membership in PCA frame · L*={args.L:g}, ab_span={args.ab_span:g}\n"
        f"white = in gamut ({100*frac_in:.1f}%)   red = out of gamut ({100*(1-frac_in):.1f}%)"
    )
    # Fixed axes rectangle so this script and render_pca_color_ramp.py
    # produce pixel-identical data areas at matching figsize/dpi.
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.10, top=0.88)

    out = pathlib.Path(args.out)
    fig.savefig(out.with_suffix(".png"), dpi=args.dpi)
    fig.savefig(out.with_suffix(".pdf"), dpi=args.dpi)
    print(f"wrote {out.with_suffix('.png')}")
    print(f"wrote {out.with_suffix('.pdf')}")

    # also report how many of the reference points are out of gamut
    A_pts = ((x1 - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    B_pts = ((x2 - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    L_pts = np.full_like(A_pts, args.L)
    lin_pts = lab_to_linear_rgb_noclip(L_pts, A_pts, B_pts)
    in_pts = np.all((lin_pts >= 0.0) & (lin_pts <= 1.0), axis=-1)
    print(f"\nGrid pixels in gamut:        {100*frac_in:.2f}%")
    print(f"Reference haplotypes in gamut: {int(in_pts.sum())}/{len(x1)} ({100*in_pts.mean():.2f}%)")


if __name__ == "__main__":
    main()
