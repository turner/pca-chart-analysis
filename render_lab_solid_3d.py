#!/usr/bin/env python3
"""
Render a single 3-D still of the sRGB color solid nested inside the
visible-color (Rösch–MacAdam optimal-color) solid, in CIELAB space.

Three things are drawn in the (a*, b*, L*) frame:

  1. "Visible-color solid"  — translucent gray outer shell: the limit of
                              human object-color perception (optimal colors
                              for a D65 illuminant, CIE 1931 2° observer).
                              CIELAB is mathematically unbounded, so THIS is
                              the meaningful "outer shell" around sRGB.
  2. "sRGB solid"           — the sRGB cube, warped into its true twisted
                              shape in Lab, faces colored by real sRGB color.
  3. "Formula region"       — the paper's color formula lives at a single
                              fixed L*, sweeping a*, b* over ±ab_span. Drawn
                              as a horizontal slice plane + bold box at that
                              L*, so you can see it reach OUTSIDE the sRGB
                              solid (the parts of the box hanging in open
                              space are colors the display cannot show).

This is the 3-D parent of the 2-D figure produced by render_lab_gamut.py:
that script draws one horizontal slice; this one shows the whole solid.

Output: a PNG (still image), no interactivity required.
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import mpl_toolkits.mplot3d.art3d as art3d
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage.color import rgb2lab
from scipy.spatial import ConvexHull

import colour


# --- CIELAB white point (D65), matching the rest of this repo --------------
_Xn, _Yn, _Zn = 0.95047, 1.00000, 1.08883
_DELTA = 6.0 / 29.0


def _f(t):
    """CIELAB forward nonlinearity."""
    return np.where(t > _DELTA ** 3, np.cbrt(t),
                    t / (3 * _DELTA ** 2) + 4.0 / 29.0)


def xyz_to_lab(XYZ):
    """XYZ (white normalized to Y=1) -> CIELAB, repo white point."""
    X, Y, Z = XYZ[..., 0], XYZ[..., 1], XYZ[..., 2]
    fx, fy, fz = _f(X / _Xn), _f(Y / _Yn), _f(Z / _Zn)
    L = 116.0 * fy - 16.0
    a = 500.0 * (fx - fy)
    b = 200.0 * (fy - fz)
    return np.stack([L, a, b], axis=-1)


# ---------------------------------------------------------------------------
# 1. Visible-color (Rösch–MacAdam optimal-color) solid
# ---------------------------------------------------------------------------
def visible_color_solid_lab(step_nm=5.0, lo=380.0, hi=720.0):
    """
    Vertices of the optimal-color solid in Lab.

    The set of all physically realizable object colors is convex in XYZ
    (it is the linear image of the convex set of reflectances R(λ) ∈ [0,1]).
    Its surface is reached by "optimal colors": reflectances that are 0/1
    with at most two transitions (a band-pass and its band-stop complement).
    We enumerate those, integrate against the D65 SPD and the CIE 1931 2°
    color-matching functions to get XYZ, then convert to Lab.
    """
    shape = colour.SpectralShape(lo, hi, step_nm)
    cmfs = colour.MSDS_CMFS["CIE 1931 2 Degree Standard Observer"].copy().align(shape)
    illum = colour.SDS_ILLUMINANTS["D65"].copy().align(shape)

    S = illum.values                       # (N,)
    xbar, ybar, zbar = cmfs.values.T       # each (N,)
    n = S.shape[0]

    # normalize so a perfect reflector (R≡1) has Y = 1 (D65 white)
    k = 1.0 / np.sum(S * ybar)
    wx, wy, wz = k * S * xbar, k * S * ybar, k * S * zbar  # per-λ weights

    # enumerate band-pass reflectances [i, j) for all i < j; complement gives
    # the band-stop family. cumulative sums make each integral O(1).
    cx = np.concatenate([[0.0], np.cumsum(wx)])
    cy = np.concatenate([[0.0], np.cumsum(wy)])
    cz = np.concatenate([[0.0], np.cumsum(wz)])
    totx, toty, totz = cx[-1], cy[-1], cz[-1]

    pts = [np.array([0.0, 0.0, 0.0]),                 # R≡0  -> black
           np.array([totx, toty, totz])]              # R≡1  -> white
    for i in range(n):
        for j in range(i + 1, n + 1):
            X = cx[j] - cx[i]
            Y = cy[j] - cy[i]
            Z = cz[j] - cz[i]
            pts.append([X, Y, Z])                      # band-pass
            pts.append([totx - X, toty - Y, totz - Z])  # band-stop complement
    XYZ = np.asarray(pts)

    hull = ConvexHull(XYZ)
    lab = xyz_to_lab(XYZ)
    # triangles as Lab vertex coords, plotted in (a, b, L) order
    tris = [lab[s][:, [1, 2, 0]] for s in hull.simplices]
    return tris


# ---------------------------------------------------------------------------
# 2. sRGB solid (warped cube) as 6 colored faces
# ---------------------------------------------------------------------------
def srgb_faces(n=24):
    """Yield (A, B, L, RGBA) for each of the 6 sRGB cube faces, in Lab."""
    u = np.linspace(0.0, 1.0, n)
    U, V = np.meshgrid(u, u)
    O = np.ones_like(U)
    Z = np.zeros_like(U)
    faces = [
        np.stack([Z, U, V], -1), np.stack([O, U, V], -1),   # R = 0, 1
        np.stack([U, Z, V], -1), np.stack([U, O, V], -1),   # G = 0, 1
        np.stack([U, V, Z], -1), np.stack([U, V, O], -1),   # B = 0, 1
    ]
    for rgb in faces:
        lab = rgb2lab(rgb[np.newaxis, ...])[0]
        rgba = np.concatenate([rgb, np.ones((*rgb.shape[:2], 1))], axis=-1)
        yield lab[..., 1], lab[..., 2], lab[..., 0], rgba


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="lab_solid_3d.png")
    ap.add_argument("--L", type=float, default=80.0,
                    help="paper's fixed L* for the formula slice")
    ap.add_argument("--ab-span", type=float, default=120.0,
                    help="formula extent: a*, b* in [-span, +span]")
    ap.add_argument("--face-res", type=int, default=24,
                    help="grid resolution per sRGB cube face")
    ap.add_argument("--elev", type=float, default=22.0)
    ap.add_argument("--azim", type=float, default=-58.0)
    ap.add_argument("--dpi", type=int, default=220)
    ap.add_argument("--no-formula", action="store_true",
                    help="omit the formula slice/box overlay")
    args = ap.parse_args()

    plt.rcParams["font.size"] = 10

    fig = plt.figure(figsize=(11, 10), dpi=args.dpi)
    ax = fig.add_subplot(111, projection="3d", computed_zorder=False)
    ax.set_facecolor("white")

    # --- outer shell: visible-color solid --------------------------------
    tris = visible_color_solid_lab()
    shell = Poly3DCollection(tris, facecolor="#7a7a7a", edgecolor="#9a9a9a",
                             linewidths=0.12, alpha=0.08)
    shell.set_zsort("min")
    ax.add_collection3d(shell)

    # --- inner solid: sRGB cube warped into Lab --------------------------
    for A, B, L, rgba in srgb_faces(args.face_res):
        ax.plot_surface(A, B, L, facecolors=rgba, rcount=args.face_res,
                        ccount=args.face_res, shade=False, antialiased=True,
                        linewidth=0)

    # --- overlay: paper's fixed-L* formula slice + box -------------------
    if not args.no_formula:
        s = args.ab_span
        L0 = args.L
        sq = np.array([[-s, -s, L0], [s, -s, L0], [s, s, L0], [-s, s, L0]])
        slab = Poly3DCollection([sq], facecolor="#ffd84a", edgecolor="none",
                                alpha=0.16)
        slab.set_zsort("max")
        ax.add_collection3d(slab)
        # bold outline of the formula extent
        loop = np.vstack([sq, sq[0]])
        ax.plot(loop[:, 0], loop[:, 1], loop[:, 2], color="#caa400",
                lw=2.2, zorder=10)

    # --- framing ---------------------------------------------------------
    lim = max(args.ab_span, 110)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(0, 100)
    ax.set_xlabel("a*  (green − / red +)", labelpad=10)
    ax.set_ylabel("b*  (blue − / yellow +)", labelpad=10)
    ax.set_zlabel("L*  (lightness)", labelpad=6)
    ax.view_init(elev=args.elev, azim=args.azim)
    ax.set_box_aspect((1, 1, 0.7))

    # legend proxies
    handles = [
        mpatches.Patch(facecolor="#7a7a7a", alpha=0.3,
                       label="Visible-color solid (human-vision limit)"),
        mpatches.Patch(facecolor="#c0506a", label="sRGB solid (displayable)"),
    ]
    if not args.no_formula:
        handles.append(mpatches.Patch(facecolor="#ffd84a", alpha=0.6,
                       label=f"Color formula slice (L*={args.L:g}, ±{args.ab_span:g})"))
    ax.legend(handles=handles, loc="upper left", framealpha=0.9, fontsize=9)

    ax.set_title(
        "The sRGB gamut is a small twisted solid inside human color vision\n"
        "CIELAB space — the paper's color formula reaches outside it",
        fontsize=12, pad=18)

    if not args.no_formula:
        fig.text(
            0.5, 0.035,
            "The gold square is the formula's color space at L*=%g. Where it "
            "extends past the\ncolored sRGB solid, those colors are clipped on "
            "screen; where it extends past the gray\nshell entirely, they are "
            "not even physically perceivable colors." % args.L,
            ha="center", va="bottom", fontsize=9, color="#333333")

    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(args.out, dpi=args.dpi, bbox_inches="tight")
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
