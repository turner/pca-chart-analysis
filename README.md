# PCLAI Color-Ramp Experiment Scripts — Cheat Sheet

Standalone scripts for exploring the failure modes described in
`pclai_color_ramp_failure_analysis.md`. All take a single TSV (with
`x1`, `x2` columns) as their only required input.

All renderers now emit **layered Photoshop (`.psd`) files** instead of
flat PNG/PDF. Each visual element lives on its own layer so it can be
toggled, re-opacity'd, or recolored in Photoshop without re-running the
script. The three PCA-frame renderers (`render_pca_psd.py`,
`render_voronoi_psd.py`, `render_kde_psd.py`) share a fixed canvas size
and axes rectangle, so their PSDs **register pixel-for-pixel** — open
them together in Photoshop and the layers stack without manual
alignment.

| Script | What it does | Output |
|---|---|---|
| **`render_pca_psd.py`** | PCA color ramp + gamut mask + reference points + axes | layered `.psd` (PCA-frame, registered) |
| `render_voronoi_psd.py` | Voronoi tessellation of superpopulation centroids over the same PCA frame | layered `.psd` (PCA-frame, registered) |
| `render_kde_psd.py` | Per-superpopulation KDE contours (95/75/50% coverage) over the same PCA frame | layered `.psd` (PCA-frame, registered) |
| `render_lab_gamut.py` | sRGB gamut at L*, reference points at TRUE a*b* + per-group Δ | layered `.psd` (CIELAB `(a*, b*)` frame) |
| `sweep_gamut.py` | Search `(L*, ab_span)` grid for in-gamut configurations | Table + recommendation to stdout |

Default parameters in the renderers reproduce the paper's published
behavior: `L*=80`, `ab_span=120`, `margin=0.1`. Override any of them.

Deprecated single-purpose flat-PNG scripts (`render_pca_color_ramp.py`,
`render_gamut_mask.py`) have been moved to `depricated/` and folded
into `render_pca_psd.py` as layers.

---

## Conda environment

A conda environment named **`pclai-chart-analysis`** is provided with
every library the scripts need. Create it once:

```bash
conda create -n pclai-chart-analysis -c conda-forge -y \
    python=3.11 numpy pandas matplotlib scipy scikit-image pip
conda run -n pclai-chart-analysis pip install pytoshop
```

(`pytoshop` is pip-only — not packaged on conda-forge.)

Then, before running any of the scripts:

```bash
conda activate pclai-chart-analysis
```

All commands below assume the env is active. Deactivate with
`conda deactivate` when done.

---

## 1. `render_pca_psd.py` — PCA color ramp + gamut mask

Output is a single Photoshop file with four registered layers
(top → bottom):

1. **Axes & labels** — matplotlib axes, ticks, PC1/PC2 labels, title with
   gamut stats.
2. **Reference points** — the haplotype scatter dots.
3. **Gamut mask** — pure red where the formula's color is out of sRGB,
   transparent where it is in gamut.
4. **Color ramp** — the paper's `lab2rgb`-clipped PCA background.

### Faithful reproduction (paper defaults)
```bash
python render_pca_psd.py --tsv reference_pca_metadata.tsv --out paper_default
```

### A fully-working configuration (100% in-gamut)
```bash
python render_pca_psd.py --tsv reference_pca_metadata.tsv \
    --L 60 --ab-span 45 --out fix_safe_L60_span45
```

### Full-opacity color ramp (no "ocean effect")
```bash
python render_pca_psd.py --tsv reference_pca_metadata.tsv \
    --bg-alpha 1.0 --out paper_default_alpha1
```

### All knobs
| flag | default | meaning |
|---|---|---|
| `--tsv` | (required) | TSV with `x1`, `x2` |
| `--out` | `pca_layered` | output basename (`.psd` appended) |
| `--L` | `80.0` | fixed CIELAB lightness |
| `--ab-span` | `120.0` | half-width of a*,b* mapping |
| `--margin` | `0.1` | padding on bounding box |
| `--bg-alpha` | `0.38` | opacity of the **Color ramp** layer (paper value) |
| `--dot-alpha` | `0.9` | opacity of the **Reference points** layer |
| `--mask-alpha` | `1.0` | opacity of the **Gamut mask** layer |
| `--dot-size` | `10.0` | scatter dot point size |
| `--grid` | `600` | background raster resolution |
| `--figwidth`/`--figheight` | `7.5`/`6.5` | canvas size in inches |
| `--dpi` | `180` | output resolution |

Stdout reports `% of grid in gamut` and `n reference haplotypes in gamut`.

PSDs are written uncompressed (~120 MB at default DPI) because
pytoshop 1.2.1's RLE codec has a `packbits` bug.

---

## 2. `render_voronoi_psd.py` — Voronoi tessellation over the PCA frame

Voronoi cells on the per-superpopulation centroids, layered on the same
canvas as `render_pca_psd.py`. Layers (top → bottom):

1. **Axes & labels**
2. **Centroids** — superpopulation centroid markers + name labels.
3. **Voronoi cells**
4. **Sample scatter** — faint per-sample dots, colored per superpopulation.

```bash
python render_voronoi_psd.py --tsv reference_pca_metadata.tsv --out voronoi_layered
```

Open alongside `pca_layered.psd` in Photoshop — every layer registers
pixel-for-pixel with the PCA chart and the KDE chart below.

### All knobs
| flag | default | meaning |
|---|---|---|
| `--tsv` | (required) | TSV with `x1`, `x2`, `Population_descriptor` |
| `--out` | `voronoi_layered` | output basename |
| `--margin` | `0.1` | padding on bounding box (must match PCA renderer) |
| `--dot-size` | `8.0` | scatter dot size |
| `--dot-alpha` | `0.4` | scatter layer opacity |
| `--line-color` | `black` | Voronoi line color |
| `--line-width` | `1.2` | Voronoi line width |
| `--line-alpha` | `0.8` | Voronoi layer opacity |
| `--centroid-size` | `120.0` | centroid marker size |
| `--figwidth`/`--figheight` | `7.5`/`6.5` | must match PCA renderer |
| `--dpi` | `180` | must match PCA renderer |

---

## 3. `render_kde_psd.py` — per-superpopulation KDE contours

Per-superpopulation density contours (95% / 75% / 50% coverage),
layered on the same canvas as `render_pca_psd.py`. Layers (top → bottom):

1. **Axes & labels**
2. **Centroids**
3. **KDE contours**
4. **Sample scatter**

```bash
python render_kde_psd.py --tsv reference_pca_metadata.tsv --out kde_layered
```

### All knobs
| flag | default | meaning |
|---|---|---|
| `--tsv` | (required) | TSV with `x1`, `x2`, `Population_descriptor` |
| `--out` | `kde_layered` | output basename |
| `--margin` | `0.1` | padding on bounding box (must match PCA renderer) |
| `--grid` | `600` | KDE evaluation grid resolution |
| `--sat`/`--val` | `0.65`/`0.85` | HSV saturation/value for contour colors |
| `--dot-size` | `4.0` | scatter dot size |
| `--dot-alpha` | `0.15` | scatter layer opacity |
| `--contour-alpha` | `1.0` | contour layer opacity |
| `--centroid-size` | `140.0` | centroid marker size |
| `--figwidth`/`--figheight` | `7.5`/`6.5` | must match PCA renderer |
| `--dpi` | `180` | must match PCA renderer |

---

## 4. `render_lab_gamut.py` — Lab a*b* plane with gamut + points

The left panel of `lab_gamut_clipping_hires.png`: sRGB gamut at the
chosen L*, the formula's ±`ab_span` box, all reference haplotypes at
their TRUE computed `(a*, b*)`, and per-group mean Δ annotations. Auto-
groups 1000G `Population_descriptor`s into AFR / EUR / EAS / SAS / AMR.

This chart is in the **CIELAB (a*, b*) frame**, not PC1/PC2, so it does
*not* register with the PCA / Voronoi / KDE PSDs. Layers (top → bottom):

1. **Axes & title**
2. **Legend**
3. **Centroids**
4. **Sample points**
5. **Origin axes**
6. **Formula bbox** — dashed yellow rectangle for ±`ab_span`.
7. **Gamut boundary** — white contour at the sRGB gamut edge.
8. **sRGB gamut** — gamut interior, colored by true sRGB; transparent outside.
9. **Background** — solid dark fill.

### Faithful reproduction
```bash
python render_lab_gamut.py --tsv reference_pca_metadata.tsv --out lab_gamut_L80
```

### See the gamut grow as L* drops
```bash
python render_lab_gamut.py --tsv reference_pca_metadata.tsv --L 70 --out lab_gamut_L70
python render_lab_gamut.py --tsv reference_pca_metadata.tsv --L 60 --out lab_gamut_L60
python render_lab_gamut.py --tsv reference_pca_metadata.tsv --L 55 --out lab_gamut_L55
```

### See points land inside the gamut when sized correctly
```bash
python render_lab_gamut.py --tsv reference_pca_metadata.tsv \
    --L 60 --ab-span 45 --out lab_gamut_fix_safe
```

### All knobs
| flag | default | meaning |
|---|---|---|
| `--tsv` | (required) | TSV with `x1`, `x2` (and ideally `Population_descriptor`) |
| `--out` | `lab_gamut` | output basename |
| `--L` | `80.0` | fixed CIELAB lightness |
| `--ab-span` | `120.0` | half-width of a*,b* mapping |
| `--margin` | `0.1` | padding on bounding box |
| `--axis-limit` | `150.0` | plot extent ± in Lab units |
| `--grid` | `600` | gamut sampling resolution |
| `--dot-size` | `4.0` | reference dot size |
| `--figwidth`/`--figheight` | `9.0`/`9.0` | figure size |
| `--dpi` | `180` | output resolution |

Prints a per-group n / meanΔ / maxΔ summary to stdout.

---

## 5. `sweep_gamut.py` — find a working `(L*, ab_span)` pair

Tests every `(L*, ab_span)` pair on a grid and reports the percentage of
reference haplotypes that fall inside the sRGB gamut for each. Prints a
heatmap-style table and recommends the configuration with the largest
`ab_span` (most saturated colors) that meets a coverage target.

```bash
python sweep_gamut.py --tsv reference_pca_metadata.tsv               # target 0.95
python sweep_gamut.py --tsv reference_pca_metadata.tsv --target 1.0  # 100% in-gamut
python sweep_gamut.py --tsv reference_pca_metadata.tsv --target 0.90
```

### All knobs
| flag | default | meaning |
|---|---|---|
| `--tsv` | (required) | TSV with `x1`, `x2` |
| `--margin` | `0.1` | padding on bounding box (must match renderers) |
| `--target` | `0.95` | required fraction of points in-gamut |

The grid (L* ∈ [40, 95] step 5, ab_span ∈ [20, 120] step 5) is hard-
coded; edit the `Ls` and `spans` lists to refine.

---

## Presentation overlay

Open `pca_layered.psd`, `voronoi_layered.psd`, and `kde_layered.psd`
together in Photoshop. Because all three were rendered against the same
`AXES_BOX`, `figsize`, and `dpi`, their layers register pixel-for-pixel.
Toggle / re-opacity any combination for the talk.

`lab_gamut_*.psd` lives in its own coordinate frame (Lab a*/b*) and is
opened separately.

---

## Empirically-verified parameter recipes

Numbers below are from `sweep_gamut.py` against
`reference_pca_metadata.tsv` (3122 haplotypes, margin=0.1).

| Goal | `--L` | `--ab-span` | % points in gamut | Notes |
|---|---|---|---|---|
| **Paper baseline** | 80 | 120 | **0.0%** | Every point clipped |
| Paper-suggested "Fix 1" | 80 | 50 | 47.3% | Improves but doesn't fix |
| Paper-suggested "Fix 2" | 55 | 90 | 28.4% | Improves but doesn't fix |
| Conservative working pick | **60** | **45** | **100.0%** | Recommended default |
| Bolder working pick | 65 | 50 | 99.3% | Slightly more saturated |
| Max saturation at L*=80 | 80 | 25 | 100% | Very washed out |
| Largest span overall | 60 | 45 | 100% | Same as conservative pick |

**Key takeaway:** at the paper's L*=80, the largest `ab_span` that keeps
100% of points in-gamut is ~25 (very pastel). To preserve saturation you
must drop to L*=60–65. The paper's joint choice of `(L*=80, ab_span=120)`
is physically incompatible with the sRGB gamut.

---

## Quick sanity check

After any `render_pca_psd.py` run, the
`Reference haplotypes in gamut: N/3122` line on stdout is the single
number that tells you whether a parameter choice actually works. Use
`sweep_gamut.py` to find working pairs systematically rather than
guessing.
