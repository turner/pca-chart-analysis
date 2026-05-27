#!/usr/bin/env python3
"""
Render the PCLAI PCA color-ramp background + reference points from a single TSV.

Reproduces the color logic of `pca_to_rgb_setup` in paintings.py but as a
standalone experiment: no results.pkl.gz, no VCFs, no PCA pickle required.

Inputs:
  --tsv          A TSV with at least columns x1, x2 (e.g. reference_pca_metadata.tsv).
                 If a Population_descriptor column is present it is used for the legend.

Knobs you'll want to sweep:
  --L            Fixed CIELAB lightness (paper uses 80).
  --ab-span      Half-width of a*,b* mapping (paper uses 120).
  --margin       Padding added to bounding box (paper uses 0.1).
  --bg-alpha     Background alpha (paper uses 0.38).
  --dot-alpha    Reference dot alpha (paper uses 0.9).
  --grid         Background grid resolution.
  --scale-by-pca PATH to pickled PCA object; if given, axes are whitened by
                 sqrt(explained_variance_) like the paper does.

Outputs:
  <out>.png and <out>.pdf
"""

import argparse
import pickle
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.color import lab2rgb


def build_color_fn(x1_min, x1_max, x2_min, x2_max, L_fixed, ab_span):
    def pca_to_rgb(x1, x2):
        a = ((x1 - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * ab_span
        b = ((x2 - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * ab_span
        L = np.full_like(a, L_fixed, dtype=float)
        Lab = np.stack([L, a, b], -1)
        return np.clip(lab2rgb(Lab[np.newaxis, ...])[0], 0, 1)
    return pca_to_rgb


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", default="pca_color_ramp")
    ap.add_argument("--L", type=float, default=80.0)
    ap.add_argument("--ab-span", type=float, default=120.0)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--bg-alpha", type=float, default=0.38)
    ap.add_argument("--dot-alpha", type=float, default=0.9)
    ap.add_argument("--dot-size", type=float, default=10.0)
    ap.add_argument("--grid", type=int, default=600)
    ap.add_argument("--scale-by-pca", default=None, help="Path to pickled PCA object (optional)")
    ap.add_argument("--figwidth", type=float, default=7.5)
    ap.add_argument("--figheight", type=float, default=6.5)
    ap.add_argument("--dpi", type=int, default=180)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    if "x1" not in df.columns or "x2" not in df.columns:
        raise SystemExit(f"TSV must contain columns x1 and x2; got {list(df.columns)}")
    x1 = df["x1"].to_numpy(float)
    x2 = df["x2"].to_numpy(float)

    if args.scale_by_pca:
        with open(args.scale_by_pca, "rb") as f:
            pca = pickle.load(f)
        s1 = float(np.sqrt(pca.explained_variance_[0]))
        s2 = float(np.sqrt(pca.explained_variance_[1]))
        x1, x2 = x1 / s1, x2 / s2

    x1_min, x1_max = x1.min() - args.margin, x1.max() + args.margin
    x2_min, x2_max = x2.min() - args.margin, x2.max() + args.margin

    pca_to_rgb = build_color_fn(x1_min, x1_max, x2_min, x2_max, args.L, args.ab_span)

    gx = np.linspace(x1_min, x1_max, args.grid)
    gy = np.linspace(x2_min, x2_max, args.grid)
    GX, GY = np.meshgrid(gx, gy)
    bg_rgb = pca_to_rgb(GX, GY)

    dot_rgb = pca_to_rgb(x1, x2)

    fig, ax = plt.subplots(figsize=(args.figwidth, args.figheight), dpi=args.dpi)
    ax.imshow(
        bg_rgb,
        extent=(x1_min, x1_max, x2_min, x2_max),
        origin="lower",
        aspect="auto",
        interpolation="bilinear",
        alpha=args.bg_alpha,
    )
    ax.scatter(
        x1, x2,
        c=dot_rgb,
        s=args.dot_size,
        alpha=args.dot_alpha,
        linewidths=0.15,
        edgecolors="black",
    )

    ax.set_xlim(x1_min, x1_max)
    ax.set_ylim(x2_min, x2_max)

    ax.set_xlabel("PC1" + (" (whitened)" if args.scale_by_pca else ""))
    ax.set_ylabel("PC2" + (" (whitened)" if args.scale_by_pca else ""))
    ax.set_title(
        f"PCLAI color ramp — L*={args.L}, ab_span={args.ab_span}, "
        f"bg_alpha={args.bg_alpha}, n={len(x1)}\n "
    )
    # Fixed axes rectangle so this script and render_gamut_mask.py
    # produce pixel-identical data areas at matching figsize/dpi.
    fig.subplots_adjust(left=0.11, right=0.97, bottom=0.10, top=0.88)

    out = pathlib.Path(args.out)
    fig.savefig(out.with_suffix(".png"), dpi=args.dpi)
    fig.savefig(out.with_suffix(".pdf"), dpi=args.dpi)
    print(f"wrote {out.with_suffix('.png')}")
    print(f"wrote {out.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
