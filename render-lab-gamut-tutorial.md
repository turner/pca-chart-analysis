# Projecting the sRGB Gamut onto the CIELAB a\*b\* Plane

This document describes, in detail, how `render_lab_gamut.py` projects the
sRGB color gamut onto a slice of CIELAB space and draws it as the colored
region in the diagram. The goal is to make the method reproducible and to
explain *why* each step is what it is — so the figure can be trusted as
evidence about the paper's color-ramp behavior.

---

## 1. What we are actually drawing

CIELAB is a 3-D color space with axes `(L*, a*, b*)`:

- `L*` — lightness, 0 (black) to 100 (white).
- `a*` — green (−) ↔ red (+).
- `b*` — blue (−) ↔ yellow (+).

The sRGB gamut — the set of colors a standard monitor can actually display
— is a *solid* sitting inside CIELAB. It is bounded because each of the
three sRGB channels is constrained to `[0, 1]`; outside that range there is
no physical pixel value that produces the color.

The paper's color formula fixes `L*` at a single value and varies only
`(a*, b*)` to encode the two PCA coordinates. So the relevant question is
**not** "what does the whole 3-D sRGB solid look like," but rather:

> At one fixed lightness `L*`, which `(a*, b*)` values correspond to colors
> sRGB can display, and which fall outside?

Geometrically, fixing `L*` slices the sRGB solid with a horizontal plane.
The intersection is a 2-D region — an irregular, roughly hexagonal patch —
in the `(a*, b*)` plane. That patch is the "sRGB gamut" the figure colors
in; everything outside it is a Lab color that the monitor cannot reproduce
and that the rendering pipeline will silently clip.

We render that slice by brute force: sample the `(a*, b*)` plane on a dense
grid, test every sample for representability in sRGB, and color it
accordingly.

> For the 3-D view — this slice shown as one cut through the full sRGB
> solid, nested inside the human-vision color solid — see
> [`render-lab-solid-3d-tutorial.md`](render-lab-solid-3d-tutorial.md) and
> `render_lab_solid_3d.py`.

---

## 2. The forward color pipeline (Lab → linear sRGB)

To test a single `(L*, a*, b*)` point we convert it all the way to *linear*
sRGB and look at whether the result stays inside `[0, 1]`. The conversion is
the standard two-stage transform CIELAB → CIE XYZ → linear sRGB.

This is implemented in `lab_to_linear_rgb_noclip` (`render_lab_gamut.py:95`):

### 2a. CIELAB → CIE XYZ

Using the standard CIELAB inverse with the D65 reference white
`(Xn, Yn, Zn) = (0.95047, 1.00000, 1.08883)`:

```
fy = (L* + 16) / 116
fx = a* / 500 + fy
fz = fy − b* / 200

X = Xn · f⁻¹(fx)
Y = Yn · f⁻¹(fy)
Z = Zn · f⁻¹(fz)
```

where the inverse nonlinearity is the usual piecewise function with
`δ = 6/29`:

```
f⁻¹(t) =  t³                         if t > δ
          3·δ²·(t − 4/29)            otherwise
```

The piecewise form (a cubic in the bulk, a linear segment near black)
matches the CIE definition exactly and avoids the numerical blow-up a bare
cube root would cause near the origin.

### 2b. CIE XYZ → linear sRGB

A single 3×3 matrix multiply with the standard sRGB/D65 primaries matrix:

```
M = [  3.24062548  −1.53720797  −0.49862860
      −0.96893071   1.87575606   0.04151752
       0.05571012  −0.20402105   1.05699594 ]

[R_lin, G_lin, B_lin]ᵀ = M · [X, Y, Z]ᵀ
```

In the code this is `XYZ @ _M.T` so it applies to a whole grid of points at
once.

### 2c. Why "noclip" matters

The function name ends in `_noclip` on purpose. The conventional display
path (`skimage`'s `lab2rgb`, used elsewhere in the script) gamma-encodes and
then **clamps** the result into `[0, 1]` so it can be shown on screen. That
clamp is exactly the failure being studied: it hides whether a color was
ever out of range.

For the *gamut test* we must look at the raw, unclamped linear RGB. A color
is representable if and only if every channel already lies in `[0, 1]`
before any clamping. Clamping first would make every point look "in range"
and destroy the measurement.

Note also that the in-gamut test is performed on **linear** RGB, before the
sRGB gamma (transfer) curve is applied. This is valid because the gamma
curve is a monotonic bijection from `[0, 1]` to `[0, 1]`: a linear value is
inside `[0, 1]` exactly when its gamma-encoded counterpart is. Testing the
linear values is simpler and avoids an unnecessary nonlinearity in the
boundary test.

---

## 3. Sampling the plane and building the mask

The projection itself is a rasterization (`render_lab_gamut.py:170`):

1. **Build a regular grid.** Over the plotted extent `±axis_limit` (default
   ±150), lay down `grid × grid` samples (default 600 × 600) in `a*` and
   `b*`:

   ```python
   ga = np.linspace(-lim, lim, grid)
   gb = np.linspace(-lim, lim, grid)
   GA, GB = np.meshgrid(ga, gb)
   GL = np.full_like(GA, L)        # constant lightness for the whole slice
   ```

   Every grid cell is one candidate Lab color `(L, GA, GB)` at the fixed
   lightness.

2. **Convert and test.** Push the entire grid through the noclip pipeline
   and mark a cell as in-gamut only if **all three** linear channels are
   within `[0, 1]`:

   ```python
   lin = lab_to_linear_rgb_noclip(GL, GA, GB)
   in_gamut = np.all((lin >= 0.0) & (lin <= 1.0), axis=-1)
   ```

   `in_gamut` is now a boolean image: `True` where this `(a*, b*)` at this
   `L*` is a real sRGB color, `False` where it is not. This boolean field
   *is* the projection of the gamut slice onto the plane.

3. **Color the interior.** Separately, run the same grid through the
   ordinary display path (`lab2rgb` + clip) to get the actual on-screen
   color of each representable cell:

   ```python
   grid_disp = np.clip(lab2rgb(stack(GL, GA, GB)), 0, 1)
   ```

   Build an RGBA raster that is fully opaque and true-colored inside the
   gamut, and fully transparent outside it (so the dark background layer
   shows through):

   ```python
   gamut_rgba = np.zeros((grid, grid, 4))
   gamut_rgba[in_gamut, :3] = grid_disp[in_gamut]
   gamut_rgba[in_gamut, 3]  = 1.0
   ```

The transparency in the alpha channel is what gives the figure its
characteristic shape: the colored blob is the gamut slice, and its silhouette
is the gamut boundary at this lightness.

> **Note — the gamut shape is emergent, never drawn.**
> There is no explicit gamut geometry anywhere in this code: no boundary
> polygon, no vertices, no parametric curve. The characteristic irregular,
> roughly hexagonal outline is the *aggregate result of ~360,000 independent
> point-inclusion tests* (one per grid cell at the default 600×600
> resolution). Each cell answers a single yes/no question — "are all three
> linear-sRGB channels within `[0, 1]`?" — and the colored blob is simply
> every cell that answered "yes," painted with its true color; every "no"
> cell is left transparent. The shape falls out of *where the answer flips*
> from yes to no across the plane; it is not constructed.
>
> This extends even to the white outline of Section 5. The boundary line is
> not an independently drawn shape either — it is a contour traced at the
> 0.5 level between the `True` and `False` cells, i.e. the edge of the same
> inclusion mask made visible. Consequently `--grid` is the *only* control
> over how finely this emergent edge is resolved: it is the sampling density
> of the inclusion test, not a smoothing or drawing parameter.

---

## 4. Placing the raster in data coordinates

The raster is drawn with `imshow` (`render_lab_gamut.py:220`):

```python
ax.imshow(gamut_rgba,
          extent=(-lim, lim, -lim, lim),
          origin="lower",
          interpolation="bilinear")
```

Two choices keep the projection geometrically honest:

- **`extent=(-lim, lim, -lim, lim)`** maps the raster's pixel grid directly
  onto the `(a*, b*)` data axes, so a sample's screen position equals its
  Lab coordinate. Combined with `ax.set_aspect("equal")`, one unit of `a*`
  equals one unit of `b*` on screen — distances and angles in the plot are
  true Lab distances.
- **`origin="lower"`** matches the orientation of `np.linspace`/`meshgrid`
  (row 0 = smallest `b*`), so the image is not flipped relative to the data.

`bilinear` interpolation only smooths the *display* of an already-computed
600×600 field; it does not change the membership test. The true boundary is
captured separately and exactly in the next step.

---

## 5. The boundary contour

The crisp white outline of the gamut is not drawn from the raster — it is
extracted from the boolean mask directly with a contour at the half-level
(`render_lab_gamut.py:232`):

```python
ax.contour(GA, GB, in_gamut.astype(float),
           levels=[0.5], colors="white", linewidths=1.5)
```

Casting the boolean mask to `0.0/1.0` and contouring at `0.5` traces exactly
the transition between in-gamut and out-of-gamut cells. This gives a clean,
resolution-faithful gamut boundary independent of the colored fill, which is
why the boundary reads sharply even where the fill fades.

---

## 6. The clipped points and their displacement paths

The figure carries each reference haplotype **twice**, on two separate,
independently toggleable layers:

- **"Sample points"** — the dot at its **true** (pre-clip) `(a*, b*)`,
  colored by its post-clip displayed sRGB.
- **"Clipped points"** — the dot at the `(a*, b*)` it **actually lands on**
  after the sRGB clip (outlined in white to distinguish it), with a third
  **"Clip paths"** layer drawing a faint connector between the two.

The landing position is computed by running each point through the ordinary
display path and back into Lab:

```python
disp_rgb    = np.clip(lab2rgb(Lab_pts), 0, 1)   # what the screen shows
landed_lab  = rgb2lab(disp_rgb)                  # that color, back in Lab
a_landed, b_landed = landed_lab[..., 1], landed_lab[..., 2]
```

### 6a. Why landing points do *not* sit on the drawn boundary

Intuition says every clipped point should snap onto the white gamut outline.
Most do not — and this is correct, not a bug. The reason is that **the clip
is a 3-D operation while the figure is a 2-D slice.**

`np.clip(rgb, 0, 1)` clamps each channel independently in RGB space. When it
pulls a channel from, say, `1.07` down to `1.0`, the resulting color differs
from the requested one in *all three* Lab dimensions — including `L*`. A
point requested at `L*=80` routinely lands at `L*=72`, or `L*=53`. It lands
honestly on the surface of the 3-D sRGB solid, but **at a different
lightness than the slice being drawn.** Projected onto the fixed-`L*=80`
plane, that landing point falls *inside* the `L*=80` boundary, because the
gamut cross-section at the lower lightness is wider and sits inboard of the
`L*=80` outline.

### 6b. The connectors are a projected shadow, not the full move

Each "Clip paths" line is the `(a*, b*)` **shadow** of a genuinely 3-D
displacement vector `(ΔL*, Δa*, Δb*)`. The `ΔL*` leg points straight out of
the page and is invisible here. Consequences worth keeping in mind when
reading the figure:

- A point that clips mostly in **lightness** shows a short on-screen
  connector despite a large true displacement — most of its move went out of
  the plane.
- The per-group **Δ** annotations (and the headline `mean Δ`) are the
  **full 3-D** Euclidean distances in Lab, so they are generally *larger*
  than the drawn connector lengths would imply.

To see the missing third leg as real depth, use the companion 3-D view
(`render_lab_solid_3d.py` / `render-lab-solid-3d-tutorial.md`), where the
clip lands visibly on the gamut surface instead of being flattened onto one
lightness plane.

---

## 7. How the projection drives the figure's claim

Everything quantitative in the figure rests on this projection:

- The reference haplotypes are placed at their **true** (pre-clip) `(a*, b*)`
  computed from the same normalization the paper uses
  (`render_lab_gamut.py:160`), so each dot's position can be compared
  directly against the gamut silhouette.
- A point is counted as out-of-gamut by indexing the very same `in_gamut`
  mask at that point's grid cell (`render_lab_gamut.py:198`) — the
  headline "`n_out`/`n_total` outside the gamut" number is read straight
  off the projection, not estimated separately.
- The dashed yellow box marks the formula's `±ab_span` extent. Seeing it
  spill far past the colored gamut blob is the visual statement of the
  failure: the formula addresses a region of Lab the display cannot show.

In short, the projection is: *fix `L*`, sample the `(a*, b*)` plane, convert
each sample to unclamped linear sRGB, and keep the cells whose channels all
land in `[0, 1]`.* That boolean field — colored inside, transparent and
white-edged at its border — is the sRGB gamut as it appears in the diagram.

---

## 8. Parameters that change the projection

| Flag | Effect on the projection |
|------|--------------------------|
| `--L` | Which horizontal slice of the sRGB solid is taken. The gamut blob's size and shape change strongly with `L*`; it is largest near mid-lightness and shrinks toward black/white. |
| `--axis-limit` | The plotted `(a*, b*)` extent (the window), `±lim`. Does not change the gamut, only how much empty space surrounds it. |
| `--grid` | Sampling resolution of the mask. Higher values sharpen the boundary contour and the in/out test at sub-unit precision; lower values are faster but blockier. |
| `--ab-span` | Not part of the gamut test — it sizes the dashed formula box and scales where the data points land, i.e. *what gets compared against* the gamut. |
| `--margin` | Expands the PCA min/max used to normalize points into `(a*, b*)`; affects point placement, not the gamut slice itself. |

The companion script `sweep_gamut.py` reuses the identical
`lab_to_linear_rgb_noclip` test across a grid of `(L*, ab_span)` values to
report the in-gamut fraction numerically, which is the tabular counterpart
to this visual projection.
