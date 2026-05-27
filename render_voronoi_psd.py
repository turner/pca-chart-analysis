#!/usr/bin/env python3
"""
Render a Voronoi-tessellation visualization as a layered .psd that registers
exactly with the PSD produced by render_pca_psd.py.

The output is a Photoshop file at the same canvas size, axes rectangle, and
data limits as render_pca_psd.py, with each visual element on its own layer
(from bottom to top):

  1. "Sample scatter"   — faint per-sample dots, coloured per superpopulation.
  2. "Voronoi cells"    — Voronoi tessellation lines on the centroids.
  3. "Centroids"        — superpopulation centroid markers + name labels.
  4. "Axes & labels"    — matplotlib axes/ticks/title rendered transparently.

All layers share one canvas at identical pixel dimensions. Open this PSD
alongside the PCA PSD (and the KDE PSD) in Photoshop and the layers will
register pixel-for-pixel.

Registration constants (must match render_pca_psd.py):
    AXES_BOX  = (0.11, 0.10, 0.86, 0.78)
    figsize   = 7.5 x 6.5
    dpi       = 180
    aspect    = "auto"
    xlim      = (x1.min - margin, x1.max + margin)
    ylim      = (x2.min - margin, x2.max + margin)
"""

import argparse
import pathlib

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.spatial import Voronoi, voronoi_plot_2d

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
    ap.add_argument("--out", default="voronoi_layered")
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--figwidth", type=float, default=7.5)
    ap.add_argument("--figheight", type=float, default=6.5)
    ap.add_argument("--dpi", type=int, default=180)
    ap.add_argument("--dot-size", type=float, default=8.0)
    ap.add_argument("--dot-alpha", type=float, default=0.4,
                    help="Opacity of the sample-scatter layer (0..1).")
    ap.add_argument("--line-color", default="black")
    ap.add_argument("--line-width", type=float, default=1.2)
    ap.add_argument("--line-alpha", type=float, default=0.8,
                    help="Opacity of the Voronoi-cells layer (0..1).")
    ap.add_argument("--centroid-size", type=float, default=120.0)
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

    centroids = df.groupby(LABEL_COL)[[X_COL, Y_COL]].mean()
    points = centroids.values
    labels = centroids.index.tolist()
    vor = Voronoi(points)

    populations = sorted(df[LABEL_COL].unique())
    cmap = plt.get_cmap("tab10")
    color_for = {sp: cmap(i) for i, sp in enumerate(populations)}

    # ---- Layer 1: sample scatter ----------------------------------------
    def draw_scatter(ax):
        for sp in populations:
            sub = df[df[LABEL_COL] == sp]
            ax.scatter(sub[X_COL], sub[Y_COL], s=args.dot_size,
                       color=color_for[sp], linewidths=0)
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    scatter_rgba = render_mpl_layer(draw_scatter, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 2: Voronoi cells -----------------------------------------
    def draw_cells(ax):
        voronoi_plot_2d(vor, ax=ax, show_points=False, show_vertices=False,
                        line_colors=args.line_color, line_width=args.line_width,
                        line_alpha=1.0)
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    cells_rgba = render_mpl_layer(draw_cells, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 3: centroids + labels ------------------------------------
    def draw_centroids(ax):
        ax.scatter(points[:, 0], points[:, 1], s=args.centroid_size, marker="X",
                   edgecolor="black", facecolor="white", linewidth=1.5, zorder=5)
        for (x, y), name in zip(points, labels):
            ax.annotate(name, (x, y), textcoords="offset points", xytext=(8, 8),
                        fontsize=11, fontweight="bold", zorder=6)
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
        ax.set_title("Voronoi tessellation of superpopulation centroids")

    axes_rgba = render_mpl_layer(draw_axes, args.figwidth, args.figheight, args.dpi)

    # ---- Assemble PSD ---------------------------------------------------
    layers = [
        rgba_to_psd_layer("Axes & labels",  axes_rgba),
        rgba_to_psd_layer("Centroids",      centroids_rgba),
        rgba_to_psd_layer("Voronoi cells",  cells_rgba,   opacity=int(round(args.line_alpha * 255))),
        rgba_to_psd_layer("Sample scatter", scatter_rgba, opacity=int(round(args.dot_alpha * 255))),
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
