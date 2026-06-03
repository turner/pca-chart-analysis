# From Failure Analysis to Working Tools — The Journey

## Starting point

This work began with the document already in this folder,
`pclai_color_ramp_failure_analysis.md`, and its companion figures. That
analysis made a strong, quantitative claim about the PCLAI paper's color
ramp: 100% of the 3,122 reference haplotypes fall outside the sRGB gamut,
the perceptual-uniformity guarantee that motivated the CIELAB choice is
voided, and the displayed color scale is geometrically a 1D projection
onto the gamut boundary instead of the intended 2D Lab surface.

The analysis was self-consistent and well-illustrated. But a written
critique — however thorough — is still a static artifact. To actually
*work with* the failure (test fixes, prepare presentations, convince
collaborators), we needed code we could run, modify, and re-render in
seconds.

## Why not just run `paintings.py`?

The paper's own `paintings.py` is not a standalone color-ramp explorer.
It's a CLI for generating chromosome paintings and PCA contour figures
from a full PCLAI inference run. To execute it we'd have needed:

- `results.pkl.gz` and `results_cp.pkl.gz` from a PCLAI run (the
  per-sample, per-chromosome, per-haplotype window predictions),
- per-chromosome VCFs to recover SNP positions in base pairs,
- a pickled fitted PCA object for the whitening step,
- and a sample ID to actually plot.

None of that is needed to study the color-ramp behavior itself. The
color logic lives in a single function — `pca_to_rgb_setup` — that takes
only a TSV of reference PCA coordinates. The rest of the program is
input plumbing and chromosome geometry, irrelevant to the failure mode.

So the first move was to extract the color logic and wrap it in a
program that runs against the one input we *do* have:
`reference_pca_metadata.tsv`.

## What we built and why

The tools accumulated in a natural sequence, each answering a question
the previous one raised.

**1. `render_pca_color_ramp.py` — reproduce the artifact.** A
faithful, parameterized reproduction of the paper's PCA color
background, with the same `lab2rgb` call and the same post-conversion
`np.clip(..., 0, 1)`. This lets us see the original chart in seconds,
and — crucially — sweep `L*`, `ab_span`, and `bg_alpha` to test the
fixes the failure analysis proposed.

**2. `render_lab_gamut.py` — make the clipping visible.** The PCA
color ramp shows colors but hides the geometry of the failure. This
script switches into the CIELAB `(a*, b*)` frame, draws the actual
sRGB gamut at the chosen lightness, the formula's ±`ab_span` bounding
box, and every reference haplotype at its *true* (pre-clip) Lab
position. A later pass added two more layers — each haplotype's
post-clip **landing** position and a connector tracing the move — which
makes the clip a visible displacement rather than just a color swap
(with the caveat that those paths are the 2-D shadow of a 3-D clip; see
the tutorial §6). It is the diagram in the failure analysis, rebuilt as
a program we can re-render with different parameters. The per-group Δ
numbers it prints to stdout match the paper's tables, confirming the
underlying math is reproduced correctly.

**3. `render_gamut_mask.py` — a clean visual for presentations.** The
Lab-frame diagram is dense. For talks and figures we wanted a single
unambiguous image: white where the formula's color is representable in
sRGB, pure red where it isn't, all in the same `(x1, x2)` frame as the
color ramp. This lets a viewer instantly understand the geometry of
what's being clipped — and, with the two scripts now layout-locked, the
mask can be overlaid directly on top of the color ramp in Photoshop or
Keynote without any registration work.

**4. `sweep_gamut.py` — turn the question into a calculation.** The
failure analysis suggested two parameter fixes (`L*=80, ab_span=50` and
`L*=55, ab_span=90`). Once we could measure in-gamut percentage in code,
the obvious next question was: do those fixes actually work? The sweep
script answers that systematically. It revealed that neither suggested
fix gets the data inside the gamut — Fix 1 lands at 47% in-gamut, Fix 2
at 28%. The real working region requires lower lightness *and* a
correspondingly modest span: `L*=60, ab_span=45` puts 100% of points
inside the gamut. The sweep turned the "what would fix it" section of
the failure analysis from intuition into a measured result.

**5. `render_pca_psd.py` — collapse the workflow into one artifact.**
With (1) and (3) we had a usable presentation recipe: render the chart,
render the mask, drop both into Photoshop as stacked layers, align by
hand (the layout-locking made the manual step trivial but it was still a
manual step). That gap — two files, one Photoshop session — was friction
every time we wanted to talk about the chart. So we combined the two
into a single program whose output *is* a layered Photoshop file. One
invocation produces a `.psd` with four registered layers — color ramp,
gamut mask, reference points, axes — stacked on a shared canvas at
pixel-identical dimensions. Toggle visibility, slide the mask's opacity,
and the failure mode reveals itself directly on top of the populations
it clips. **This is now the preferred renderer.** `render_pca_color_ramp.py`
and `render_gamut_mask.py` remain in the tree for archival reference but
are deprecated; new work should go through `render_pca_psd.py`.

The conceptual move at this step was small but useful: the chart and
its failure projection are not two separate visualizations that happen
to be informative when overlaid — they're one visualization, and the
artifact that carries the discussion should reflect that. A `.psd` with
both layers is harder to show in a paper, but easier to show in a talk,
in a Slack thread, or in a working session with collaborators — which
is where most of the actual conversation about the failure happens.

## A couple of representative uses

**Reproduce the paper and inspect a candidate fix side by side, as
layered Photoshop files:**
```bash
python render_pca_psd.py --tsv reference_pca_metadata.tsv \
    --out paper_default
python render_pca_psd.py --tsv reference_pca_metadata.tsv \
    --L 60 --ab-span 45 --out fixed
```
Each `.psd` carries its own chart and its own gamut projection on
registered layers — open one, hide everything except `Color ramp` for a
clean chart, then re-enable `Gamut mask` to show exactly which
populations are being clipped.

**Make a presentation slide where the gamut-mask "red" hovers over the
color ramp at exactly the right populations:** the layered `.psd` already
does this. Open it in Photoshop and the `Gamut mask` layer is registered
pixel-perfectly on top of `Color ramp` — no manual overlay step.

Operational detail is in `README.md`.

## The point of all of this

The failure analysis identified a real, quantifiable defect in a
published method. But for that finding to be useful — to share with
collaborators, present at a meeting, propose a concrete fix to the
PCLAI authors, or rebuild the visualization in a downstream paper —
the analysis needs to be *executable*. These scripts (now centered on
`render_pca_psd.py`, with `sweep_gamut.py` and `render_lab_gamut.py`
as supporting tools) turn the document's claims into things you can
run, see, modify, and verify.

Specifically:

- **Trust** — the per-group Δ numbers and the 16.4% gamut-area figure
  come out of the code matching the failure analysis to one decimal,
  confirming the analysis is computationally correct.
- **Exploration** — every parameter the analysis discusses (L*,
  `ab_span`, `margin`, `bg_alpha`) is now a flag. Hypotheses about
  what would help can be tested in seconds, not argued about in
  abstract.
- **Communication** — the layered `.psd` and the sweep table give you,
  for a presentation, exactly the artifacts an audience needs: one
  file whose layers can be toggled to show what the paper produces,
  what's structurally wrong with it, and (re-run with a working
  `L*`/`ab_span`) what a corrected version looks like.

The analysis told us *what* was broken. These tools let us show *how*
it's broken, *how much*, and *what to do about it*.
