#!/usr/bin/env python3
"""
Render the LEFT panel of the paper's gamut-clipping figure as a layered .psd.

Plots, in the CIELAB a*b* plane at fixed L*, with each visual element on its
own layer (bottom to top):

  1. "Background"        — solid dark fill behind everything.
  2. "sRGB gamut"        — sRGB gamut interior, colored by true sRGB;
                           transparent outside the gamut.
  3. "Gamut boundary"    — white contour at the sRGB gamut edge.
  4. "Formula bbox"      — dashed yellow rectangle for the formula extent
                           a*, b* in [-ab_span, +ab_span].
  5. "Origin axes"       — faint white lines through (0, 0).
  6. "Sample points"     — reference haplotypes at their TRUE (a*, b*),
                           colored by their POST-CLIP displayed sRGB.
  7. "Centroids"         — per-group centroid markers + name/Δ labels.
  8. "Legend"            — corner annotation box.
  9. "Axes & title"      — matplotlib axes/ticks/labels/title.

Input: a TSV with x1, x2 (e.g. reference_pca_metadata.tsv).
If a Population_descriptor column is present, points are grouped into
continental groups using a built-in 1000G mapping; otherwise everything
is one group.

The (a*, b*) coordinates are computed with the same normalization the
paper uses:
  a* = ((x1 - x1_min)/(x1_max - x1_min) - 0.5) * 2 * ab_span
  b* = ((x2 - x2_min)/(x2_max - x2_min) - 0.5) * 2 * ab_span
with x1_min/max expanded by --margin.
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from skimage.color import lab2rgb, rgb2lab

from pytoshop.user import nested_layers
from pytoshop import enums


# --- 1000G population descriptor -> continental group ---------------------
GROUP_MAP = {
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

GROUP_COLORS = {
    "AFR": "#f0a020",
    "EUR": "#e83e8c",
    "EAS": "#17a2b8",
    "SAS": "#9b59b6",
    "AMR": "#27ae60",
}

BG_COLOR = "#0a0a14"
AXES_BOX = (0.08, 0.08, 0.88, 0.84)


# --- Lab -> linear sRGB, NO clipping (for gamut testing) ------------------
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
    ap.add_argument("--out", default="lab_gamut")
    ap.add_argument("--L", type=float, default=80.0)
    ap.add_argument("--ab-span", type=float, default=120.0)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--axis-limit", type=float, default=150.0, help="a*,b* plot extent (±)")
    ap.add_argument("--grid", type=int, default=600, help="resolution of gamut sampling grid")
    ap.add_argument("--dot-size", type=float, default=4.0)
    ap.add_argument("--figwidth", type=float, default=9.0)
    ap.add_argument("--figheight", type=float, default=9.0)
    ap.add_argument("--dpi", type=int, default=180)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    if "x1" not in df.columns or "x2" not in df.columns:
        raise SystemExit(f"TSV must contain columns x1 and x2; got {list(df.columns)}")

    x1 = df["x1"].to_numpy(float)
    x2 = df["x2"].to_numpy(float)

    x1_min, x1_max = x1.min() - args.margin, x1.max() + args.margin
    x2_min, x2_max = x2.min() - args.margin, x2.max() + args.margin

    a_true = ((x1 - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    b_true = ((x2 - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    L_arr = np.full_like(a_true, args.L)

    Lab_pts = np.stack([L_arr, a_true, b_true], -1)[np.newaxis, ...]
    disp_rgb = np.clip(lab2rgb(Lab_pts)[0], 0, 1)

    landed_lab = rgb2lab(disp_rgb[np.newaxis, ...])[0]
    delta = np.sqrt(np.sum((landed_lab - np.stack([L_arr, a_true, b_true], -1)) ** 2, axis=-1))

    # --- Build gamut raster at this L* -----------------------------------
    lim = args.axis_limit
    ga = np.linspace(-lim, lim, args.grid)
    gb = np.linspace(-lim, lim, args.grid)
    GA, GB = np.meshgrid(ga, gb)
    GL = np.full_like(GA, args.L)

    lin = lab_to_linear_rgb_noclip(GL, GA, GB)
    in_gamut = np.all((lin >= 0.0) & (lin <= 1.0), axis=-1)

    grid_disp = lab2rgb(np.stack([GL, GA, GB], -1)[np.newaxis, ...])[0]
    grid_disp = np.clip(grid_disp, 0, 1)

    # RGBA canvas: transparent outside the gamut so the Background layer shows through
    gamut_rgba = np.zeros((args.grid, args.grid, 4), dtype=float)
    gamut_rgba[in_gamut, :3] = grid_disp[in_gamut]
    gamut_rgba[in_gamut, 3] = 1.0

    # group assignment
    if "Population_descriptor" in df.columns:
        groups = df["Population_descriptor"].map(GROUP_MAP).fillna("OTHER").to_numpy()
    else:
        groups = np.full(len(df), "ALL", dtype=object)

    unique_groups = [g for g in ["AFR", "EUR", "EAS", "SAS", "AMR"] if (groups == g).any()]
    other_present = (groups == "OTHER").any() or (groups == "ALL").any()

    n_total = len(df)
    n_out = int((~in_gamut[
        np.clip(((b_true + lim) / (2 * lim) * (args.grid - 1)).astype(int), 0, args.grid - 1),
        np.clip(((a_true + lim) / (2 * lim) * (args.grid - 1)).astype(int), 0, args.grid - 1),
    ]).sum())

    def setup_axes(ax):
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.set_axis_off()

    # ---- Layer 1: solid background ---------------------------------------
    def draw_background(ax):
        ax.add_patch(mpatches.Rectangle(
            (-lim, -lim), 2 * lim, 2 * lim,
            facecolor=BG_COLOR, edgecolor="none", zorder=0,
        ))
        setup_axes(ax)

    bg_rgba = render_mpl_layer(draw_background, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 2: sRGB gamut raster (transparent outside) ---------------
    def draw_gamut(ax):
        ax.imshow(
            gamut_rgba,
            extent=(-lim, lim, -lim, lim),
            origin="lower",
            interpolation="bilinear",
        )
        setup_axes(ax)

    gamut_layer_rgba = render_mpl_layer(draw_gamut, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 3: gamut boundary contour --------------------------------
    def draw_boundary(ax):
        ax.contour(
            GA, GB, in_gamut.astype(float),
            levels=[0.5], colors="white", linewidths=1.5,
        )
        setup_axes(ax)

    boundary_rgba = render_mpl_layer(draw_boundary, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 4: formula bounding box ----------------------------------
    s = args.ab_span
    def draw_bbox(ax):
        ax.add_patch(mpatches.Rectangle(
            (-s, -s), 2 * s, 2 * s,
            fill=False, edgecolor="#ffd84a", linestyle="--", linewidth=1.4,
        ))
        setup_axes(ax)

    bbox_rgba = render_mpl_layer(draw_bbox, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 5: origin axes -------------------------------------------
    def draw_origin_axes(ax):
        ax.axhline(0, color="white", linewidth=0.4, alpha=0.5)
        ax.axvline(0, color="white", linewidth=0.4, alpha=0.5)
        setup_axes(ax)

    origin_rgba = render_mpl_layer(draw_origin_axes, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 6: sample points -----------------------------------------
    def draw_points(ax):
        for g in unique_groups:
            mask = groups == g
            ax.scatter(
                a_true[mask], b_true[mask],
                c=disp_rgb[mask],
                s=args.dot_size,
                linewidths=0.0,
                alpha=0.85,
            )
        if other_present:
            mask = (groups == "OTHER") | (groups == "ALL")
            ax.scatter(
                a_true[mask], b_true[mask],
                c=disp_rgb[mask], s=args.dot_size,
                linewidths=0.0, alpha=0.85,
            )
        setup_axes(ax)

    points_rgba = render_mpl_layer(draw_points, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 7: centroids + group labels ------------------------------
    def draw_centroids(ax):
        for g in unique_groups:
            mask = groups == g
            ca, cb = a_true[mask].mean(), b_true[mask].mean()
            d = delta[mask].mean()
            ax.scatter([ca], [cb], s=120, facecolor="none",
                       edgecolor=GROUP_COLORS[g], linewidths=2.0, zorder=5)
            ax.annotate(
                f"{g}\n(Δ={d:.0f})",
                xy=(ca, cb), xytext=(12, 12), textcoords="offset points",
                color=GROUP_COLORS[g], fontsize=11, fontweight="bold",
                ha="left", va="bottom",
            )
        setup_axes(ax)

    centroids_rgba = render_mpl_layer(draw_centroids, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 8: legend annotation -------------------------------------
    def draw_legend(ax):
        ax.text(
            -lim + 6, -lim + 6,
            "white contour = sRGB gamut boundary at this L*\n"
            f"dashed yellow box = formula extent ±{s:g}\n"
            "dots at TRUE (a*, b*); colored by POST-CLIP sRGB",
            color="#ffe680", fontsize=8.5,
            va="bottom", ha="left",
            bbox=dict(boxstyle="round,pad=0.4", facecolor=BG_COLOR,
                      edgecolor="#ffd84a", alpha=0.85),
        )
        setup_axes(ax)

    legend_rgba = render_mpl_layer(draw_legend, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 9: axes & title (visible spines/ticks/labels) ------------
    title = (
        f"PCLAI Color Formula: sRGB Gamut Clipping at L* = {args.L:g}\n"
        f"{n_out}/{n_total} ({100*n_out/n_total:.1f}%) reference haplotypes fall outside the sRGB gamut · "
        f"mean Δ = {delta.mean():.1f} Lab units"
    )

    def draw_axes(ax):
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_aspect("equal")
        ax.tick_params(colors="white", labelsize=9)
        for spine in ax.spines.values():
            spine.set_color("white")
        ax.set_xlabel("a*   (green ← 0 → red)", color="white", fontsize=11)
        ax.set_ylabel("b*   (blue ← 0 → yellow)", color="white", fontsize=11)
        ax.set_title(title, color="white", fontsize=11, pad=14)

    axes_rgba = render_mpl_layer(draw_axes, args.figwidth, args.figheight, args.dpi)

    # ---- Assemble PSD (top -> bottom in pytoshop list order) -----------
    layers = [
        rgba_to_psd_layer("Axes & title",   axes_rgba),
        rgba_to_psd_layer("Legend",         legend_rgba),
        rgba_to_psd_layer("Centroids",      centroids_rgba),
        rgba_to_psd_layer("Sample points",  points_rgba),
        rgba_to_psd_layer("Origin axes",    origin_rgba),
        rgba_to_psd_layer("Formula bbox",   bbox_rgba),
        rgba_to_psd_layer("Gamut boundary", boundary_rgba),
        rgba_to_psd_layer("sRGB gamut",     gamut_layer_rgba),
        rgba_to_psd_layer("Background",     bg_rgba),
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
    print(f"  canvas: {bg_rgba.shape[1]} x {bg_rgba.shape[0]} px")

    # summary
    print(f"\nSummary at L*={args.L}, ab_span={args.ab_span}:")
    print(f"  total haplotypes: {n_total}")
    print(f"  out of gamut:     {n_out} ({100*n_out/n_total:.1f}%)")
    print(f"  a* range:         [{a_true.min():.1f}, {a_true.max():.1f}]")
    print(f"  b* range:         [{b_true.min():.1f}, {b_true.max():.1f}]")
    print(f"  mean Δ:           {delta.mean():.1f}  median: {np.median(delta):.1f}  max: {delta.max():.1f}")
    for g in unique_groups:
        m = groups == g
        print(f"  {g}: n={m.sum():4d}  meanΔ={delta[m].mean():5.1f}  maxΔ={delta[m].max():5.1f}")


if __name__ == "__main__":
    main()
