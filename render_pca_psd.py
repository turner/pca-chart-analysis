#!/usr/bin/env python3
"""
Render the PCLAI PCA color ramp + sRGB gamut mask as a single layered .psd.

Combines the logic of `render_pca_color_ramp.py` and `render_gamut_mask.py`
into one program whose output is a Photoshop file with each visualization on
its own layer (from bottom to top):

  1. "Color ramp"       — the paper's PCA color background (lab2rgb + clip).
  2. "Gamut mask"       — pure red where the formula's color is out of sRGB,
                          transparent where it's in gamut.
  3. "Reference points" — the haplotype scatter dots, coloured per the paper.
  4. "Axes & labels"    — matplotlib axes/ticks/title rendered transparently.

All four layers share one canvas at identical pixel dimensions, with the data
area locked to the same axes rectangle used by the existing PNG tools, so the
layers register exactly on top of one another. In Photoshop you can toggle
visibility, change layer opacities, etc.

Inputs match the two source scripts:
  --tsv         TSV with x1, x2 (e.g. reference_pca_metadata.tsv)
  --L           Fixed CIELAB lightness (paper uses 80).
  --ab-span     Half-width of a*, b* mapping (paper uses 120).
  --margin      Padding around bounding box (paper uses 0.1).
  --bg-alpha    Layer opacity for the color ramp layer (paper uses 0.38).
  --dot-alpha   Layer opacity for the reference-point layer (paper uses 0.9).
  --grid        Background grid resolution.

Output:
  <out>.psd    Layered Photoshop file.
"""

import argparse
import pathlib

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from skimage.color import lab2rgb

from pytoshop.user import nested_layers
from pytoshop import enums


# ---- Lab -> linear sRGB without clipping (for the gamut mask) -------------
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


# ---- Helpers --------------------------------------------------------------

# Match the axes rectangle used by render_pca_color_ramp.py and
# render_gamut_mask.py: left=0.11, right=0.97, bottom=0.10, top=0.88.
AXES_BOX = (0.11, 0.10, 0.97 - 0.11, 0.88 - 0.10)  # left, bottom, w, h


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


# ---- Main -----------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--out", default="pca_layered")
    ap.add_argument("--L", type=float, default=80.0)
    ap.add_argument("--ab-span", type=float, default=120.0)
    ap.add_argument("--margin", type=float, default=0.1)
    ap.add_argument("--bg-alpha", type=float, default=0.38,
                    help="Opacity of the color-ramp layer (0..1).")
    ap.add_argument("--dot-alpha", type=float, default=0.9,
                    help="Opacity of the reference-points layer (0..1).")
    ap.add_argument("--mask-alpha", type=float, default=1.0,
                    help="Opacity of the gamut-mask layer (0..1).")
    ap.add_argument("--dot-size", type=float, default=10.0)
    ap.add_argument("--grid", type=int, default=600)
    ap.add_argument("--figwidth", type=float, default=7.5)
    ap.add_argument("--figheight", type=float, default=6.5)
    ap.add_argument("--dpi", type=int, default=180)
    args = ap.parse_args()

    df = pd.read_csv(args.tsv, sep="\t")
    if "x1" not in df.columns or "x2" not in df.columns:
        raise SystemExit(f"TSV must contain columns x1 and x2; got {list(df.columns)}")
    x1 = df["x1"].to_numpy(float)
    x2 = df["x2"].to_numpy(float)

    x1_min, x1_max = x1.min() - args.margin, x1.max() + args.margin
    x2_min, x2_max = x2.min() - args.margin, x2.max() + args.margin

    # ---- Layer 1: color ramp background ----------------------------------
    def pca_to_rgb(xx1, xx2):
        a = ((xx1 - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
        b = ((xx2 - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
        L = np.full_like(a, args.L, dtype=float)
        Lab = np.stack([L, a, b], -1)
        return np.clip(lab2rgb(Lab[np.newaxis, ...])[0], 0, 1)

    gx = np.linspace(x1_min, x1_max, args.grid)
    gy = np.linspace(x2_min, x2_max, args.grid)
    GX, GY = np.meshgrid(gx, gy)
    bg_rgb = pca_to_rgb(GX, GY)
    dot_rgb = pca_to_rgb(x1, x2)

    def draw_ramp(ax):
        ax.imshow(
            bg_rgb,
            extent=(x1_min, x1_max, x2_min, x2_max),
            origin="lower",
            aspect="auto",
            interpolation="bilinear",
        )
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    ramp_rgba = render_mpl_layer(draw_ramp, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 2: gamut mask --------------------------------------------
    A = ((GX - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    B = ((GY - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    Lfield = np.full_like(A, args.L)
    lin = lab_to_linear_rgb_noclip(Lfield, A, B)
    in_gamut = np.all((lin >= 0.0) & (lin <= 1.0), axis=-1)
    frac_in = float(in_gamut.mean())

    mask_img = np.zeros((args.grid, args.grid, 4), dtype=float)
    mask_img[~in_gamut] = [1.0, 0.0, 0.0, 1.0]   # opaque red where out of gamut
    # in-gamut pixels stay (0,0,0,0) -> transparent

    def draw_mask(ax):
        ax.imshow(
            mask_img,
            extent=(x1_min, x1_max, x2_min, x2_max),
            origin="lower",
            aspect="auto",
            interpolation="nearest",
        )
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    mask_rgba = render_mpl_layer(draw_mask, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 3: reference points --------------------------------------
    def draw_points(ax):
        ax.scatter(
            x1, x2,
            c=dot_rgb,
            s=args.dot_size,
            linewidths=0.15,
            edgecolors="black",
        )
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_axis_off()

    pts_rgba = render_mpl_layer(draw_points, args.figwidth, args.figheight, args.dpi)

    # ---- Layer 4: axes & labels -----------------------------------------
    A_pts = ((x1 - x1_min) / (x1_max - x1_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    B_pts = ((x2 - x2_min) / (x2_max - x2_min + 1e-12) - 0.5) * 2.0 * args.ab_span
    L_pts = np.full_like(A_pts, args.L)
    lin_pts = lab_to_linear_rgb_noclip(L_pts, A_pts, B_pts)
    in_pts = np.all((lin_pts >= 0.0) & (lin_pts <= 1.0), axis=-1)

    def draw_axes(ax):
        ax.set_xlim(x1_min, x1_max)
        ax.set_ylim(x2_min, x2_max)
        ax.set_xlabel("PC1")
        ax.set_ylabel("PC2")
        ax.set_title(
            f"PCLAI color ramp + sRGB gamut mask · L*={args.L:g}, "
            f"ab_span={args.ab_span:g}\n"
            f"grid in gamut: {100*frac_in:.1f}%   "
            f"haplotypes in gamut: {int(in_pts.sum())}/{len(x1)} "
            f"({100*in_pts.mean():.1f}%)"
        )

    axes_rgba = render_mpl_layer(draw_axes, args.figwidth, args.figheight, args.dpi)

    # ---- Assemble PSD ---------------------------------------------------
    layers = [
        rgba_to_psd_layer("Axes & labels",    axes_rgba),
        rgba_to_psd_layer("Reference points", pts_rgba,  opacity=int(round(args.dot_alpha * 255))),
        rgba_to_psd_layer("Gamut mask",       mask_rgba, opacity=int(round(args.mask_alpha * 255))),
        rgba_to_psd_layer("Color ramp",       ramp_rgba, opacity=int(round(args.bg_alpha * 255))),
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
    print(f"  canvas: {ramp_rgba.shape[1]} x {ramp_rgba.shape[0]} px")
    print(f"  grid pixels in gamut:         {100*frac_in:.2f}%")
    print(f"  reference haplotypes in gamut: {int(in_pts.sum())}/{len(x1)} "
          f"({100*in_pts.mean():.2f}%)")


if __name__ == "__main__":
    main()
