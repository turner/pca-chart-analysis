#!/usr/bin/env python3
"""
Render per-superpopulation KDE contours as a layered .psd that registers
exactly with the PSD produced by render_pca_psd.py.

The output is a Photoshop file at the same canvas size, axes rectangle, and
data limits as render_pca_psd.py, with each visual element on its own layer
(from bottom to top):

  1. "Sample scatter"   — faint per-sample dots, coloured per superpopulation.
  2. "KDE contours"     — 95% / 75% / 50% coverage contours per superpop.
  3. "Centroids"        — superpopulation centroid markers + name labels.
  4. "Axes & labels"    — matplotlib axes/ticks/title rendered transparently.

All layers share one canvas at identical pixel dimensions. Open this PSD
alongside the PCA PSD (and the Voronoi PSD) in Photoshop and the layers
will register pixel-for-pixel.

Registration constants (must match render_pca_psd.py):
    AXES_BOX  = (0.11, 0.10, 0.86, 0.78)
    figsize   = 7.5 x 6.5
    dpi       = 180
    aspect    = "auto"
    xlim      = (x1.min - margin, x1.max + margin)
    ylim      = (x2.min - margin, x2.max + margin)
"""

import argparse
import colorsys
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from pytoshop.user import nested_layers
from pytoshop import enums


POP_COL = "Population_descriptor"
X_COL, Y_COL = "x1", "x2"
LABEL_COL = "superpopulation"

POP_TO_SUPERPOP = {
    "Esan in Nigeria": "AFR",
    "Gambian in Western Division Ð Mandinka": "AFR",
    "Luhya in Webuye, Kenya": "AFR",
    "Mende in Sierra Leone": "AFR",
    "Yoruba in Ibadan, Nigeria": "AFR",
    "British from England and Scotland": "EUR",
    "Finnish in Finland": "EUR",
    "Iberian Populations in Spain": "EUR",
    "Toscani in Italia": "EUR",
    "Chinese Dai in Xishuangbanna": "EAS",
    "Han Chinese in Beijing, China": "EAS",
    "Han Chinese South, China": "EAS",
    "Japanese in Tokyo, Japan": "EAS",
    "Kinh in Ho Chi Minh City, Vietnam": "EAS",
    "Bengali in Bangladesh": "SAS",
    "Gujarati Indians in Houston, Texas, USA": "SAS",
    "Indian Telugu in the UK": "SAS",
    "Punjabi in Lahore, Pakistan": "SAS",
    "Sri Lankan Tamil in the UK": "SAS",
    "Mexican Ancestry in Los Angeles, California, USA": "AMR",
    "Peruvian in Lima, Peru": "AMR",
}

# These must remain identical to render_pca_psd.py.
AXES_BOX = (0.11, 0.10, 0.97 - 0.11, 0.88 - 0.10)

COVERAGE_LEVELS = [0.95, 0.75, 0.50]


def render_mpl_layer(draw_fn, figwidth, figheight, dpi):
    """Run draw_fn(ax) on a transparent figure with the shared axes box,
    return an (H, W, 4) uint8 RGBA array."""
    fig = plt.figure(figsize=(figwidth, figheight), dpi=dpi)
    fig.patch.set_alpha(0.0)
    ax = fig.add_axes(AXES_BOX)
    ax.patch.set_alpha(0.0)
    draw_fn(ax)
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba()).copy()
    plt.close(fig)
    return buf


def rgba_to_psd_layer(name, rgba, opacity=255):
    """Wrap an (H, W, 4) uint8 RGBA array as a pytoshop nested_layers.Image."""
    r = np.ascontiguousarray(rgba[..., 0])
    g = np.ascontiguousarray(rgba[..., 1])
    b = np.ascontiguousarray(rgba[..., 2])
    a = np.ascontiguousarray(rgba[..., 3])
    return nested_layers.Image(
        name=name,
        visible=True,
        opacity=int(opacity),
        color_mode=enums.ColorMode.rgb,
        channels={0: r, 1: g, 2: b, -1: a},
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", default="kde_layered")
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--figwidth", type=float, default=7.5)
    ap.add_argument("--figheight", type=float, default=6.5)
    ap.add_argument("--dpi", type=int, default=180)
    ap.add_argument("--grid", type=int, default=600)
    ap.add_argument("--sat", type=float, default=0.65)
    ap.add_argument("--val", type=float, default=0.85)
    ap.add_argument("--dot-size", type=float, default=4.0)
    ap.add_argument("--dot-alpha", type=float, default=0.15,
                    help="Opacity of the sample-scatter layer (0..1).")
    ap.add_argument("--contour-alpha", type=float, default=1.0,
                    help="Opacity of the KDE-contours layer (0..1).")
    ap.add_argument("--centroid-size", type=float, default=140.0)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    df[LABEL_COL] = df[POP_COL].map(POP_TO_SUPERPOP)
    missing = df[df[LABEL_COL].isna()][POP_COL].unique()
    if len(missing):
        raise SystemExit(f"Unmapped populations: {list(missing)}")

    x1 = df[X_COL].to_numpy(float)
    x2 = df[Y_COL].to_numpy(float)
    x1_min, x1_max = x1.min() - args.margin, x1.max() + args.margin
    x2_min, x2_max = x2.min() - args.margin, x2.max() + args.margin

    populations = sorted(df[LABEL_COL].unique())
    n = len(populations)
    palette = {sp: colorsys.hsv_to_rgb(i / n, args.sat, args.val)
               for i, sp in enumerate(populations)}

    xs = np.linspace(x1_min, x1_max, args.grid)
    ys = np.linspace(x2_min, x2_max, args.grid)
    XX, YY = np.meshgrid(xs, ys)
    grid = np.vstack([XX.ravel(), YY.ravel()])

    centroids = df.groupby(LABEL_COL)[[X_COL, Y_COL]].mean()

    # ---- Layer 1: sample scatter ----------------------------------------
    def draw_scatter(ax):
        for sp in populations:
            sub = df[df[LABEL_COL] == sp]
            ax.scatter(sub[X_COL], sub[Y_COL], s=args.dot_size,
                       color=palette[sp], linewidths=0)
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    scatter_rgba = render_mpl_layer(draw_scatter, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 2: KDE contours ------------------------------------------
    def draw_contours(ax):
        for sp in populations:
            sub = df[df[LABEL_COL] == sp][[X_COL, Y_COL]].values
            if len(sub) < 5:
                continue
            try:
                kde = gaussian_kde(sub.T)
            except (np.linalg.LinAlgError, ValueError):
                continue
            sample_density = kde(sub.T)
            sorted_d = np.sort(sample_density)[::-1]
            cdf = np.cumsum(sorted_d) / sorted_d.sum()
            thresholds = [sorted_d[min(np.searchsorted(cdf, c), len(sorted_d) - 1)]
                          for c in COVERAGE_LEVELS]
            thresholds = sorted(thresholds)
            Z = kde(grid).reshape(XX.shape)

            color = palette[sp]
            t95, t75, t50 = thresholds
            ax.contour(XX, YY, Z, levels=[t95], colors=[color],
                       linewidths=1.8, linestyles="solid")
            ax.contour(XX, YY, Z, levels=[t75], colors=[color],
                       linewidths=1.3, linestyles="dashed")
            ax.contour(XX, YY, Z, levels=[t50], colors=[color],
                       linewidths=1.0, linestyles="dotted")
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    contours_rgba = render_mpl_layer(draw_contours, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 3: centroids + labels ------------------------------------
    def draw_centroids(ax):
        for sp, (cen_x, cen_y) in centroids.iterrows():
            ax.scatter([cen_x], [cen_y], s=args.centroid_size, marker="X",
                       facecolor="white", edgecolor=palette[sp],
                       linewidths=2.0, zorder=5)
            ax.annotate(sp, (cen_x, cen_y), xytext=(8, 8),
                        textcoords="offset points",
                        fontsize=12, fontweight="bold",
                        color=tuple(c * 0.6 for c in palette[sp]), zorder=6)
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    centroids_rgba = render_mpl_layer(draw_centroids, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 4: axes & labels -----------------------------------------
    def draw_axes(ax):
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title("Per-superpopulation KDE contours (95% / 75% / 50%)")

    axes_rgba = render_mpl_layer(draw_axes, args.figwidth, args.figheight, args.dpi)

    # ---- Assemble PSD ---------------------------------------------------
    layers = [
        rgba_to_psd_layer("Axes & labels",  axes_rgba),
        rgba_to_psd_layer("Centroids",      centroids_rgba),
        rgba_to_psd_layer("KDE contours",   contours_rgba, opacity=int(round(args.contour_alpha * 255))),
        rgba_to_psd_layer("Sample scatter", scatter_rgba,  opacity=int(round(args.dot_alpha * 255))),
    ]

    psd = nested_layers.nested_layers_to_psd(
        layers,
        color_mode=enums.ColorMode.rgb,
        compression=enums.Compression.raw,
    )

    out = pathlib.Path(args.out).with_suffix(".psd")
    with open(out, "wb") as f:
        psd.write(f)

    print(f"wrote {out}")
    print(f"  canvas: {scatter_rgba.shape[1]} x {scatter_rgba.shape[0]} px")


if __name__ == "__main__":
    main()
