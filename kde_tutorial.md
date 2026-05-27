# KDE Outline Contours: A Tutorial

This document explains the algorithm used in `points2kde.py`. The goal is to render each superpopulation as a set of nested outline rings on the (x1, x2) plane — rings whose interpretation is "X% of this population's probability mass lies inside this curve."

The technique has three pieces:

1. **Kernel density estimation (KDE)** — turn a finite set of sample points into a smooth probability density.
2. **Coverage-level thresholding** — pick the density values whose level sets enclose a target fraction of the mass (50%, 75%, 95%).
3. **Contour rendering** — evaluate the density on a grid and draw the level sets as curves.

Each piece is independent and worth understanding on its own.

---

## 1. Kernel density estimation

### The problem

You have N sample points in 2D — say, the (x1, x2) coordinates of every individual labeled `EUR`. You want a smooth, continuous estimate of the probability density `p(x, y)` from which those points were drawn. With that density in hand you can:

- ask "how dense is this population at point (x, y)?"
- draw level sets — contours where the density takes a specified value
- compute coverage regions ("the smallest region containing 95% of the mass")

### The kernel idea

KDE answers this with a simple construction: place a small bump (the "kernel") at every sample point, then sum the bumps and renormalize. The result is a smooth density that is high where samples cluster and low where they are sparse.

```
density at point x  =  (1/N) * sum over i of  K_h( x - x_i )
```

where `x_i` are the sample points and `K_h` is a kernel function with bandwidth `h`. With a Gaussian kernel, `K_h` is a 2D Gaussian centered at zero with covariance scaled by `h`.

Mental picture: drop a tiny Gaussian hill on top of each data point. Add up all the hills. Where points are dense, hills overlap and stack into a ridge. Where points are isolated, you get a single small bump. The result is a smooth surface over the plane.

### Bandwidth

The bandwidth `h` controls how wide each kernel is — and therefore how smooth the resulting estimate is.

- **Small h**: each bump is narrow. The density follows the samples closely and looks lumpy. You can see individual points.
- **Large h**: each bump is wide. The density is smooth but blurs out real structure.

Choosing `h` is the only real hyperparameter in this method. `scipy.stats.gaussian_kde` uses **Scott's rule** by default — `h ∝ N^(-1/(d+4))`, where d is the dimensionality. For 2D and a few hundred points, Scott's rule produces a reasonable middle-of-the-road bandwidth without manual tuning. (The KDE memo notes this as an open question — for elongated populations a per-population or covariance-aware bandwidth would be more faithful.)

### What you have after step 1

A function `kde(x, y)` that returns a density value at any point in the plane. It's a Python callable; under the hood it's a sum of Gaussians, one per sample. You can evaluate it at the original sample points, on a regular grid, or anywhere else.

---

## 2. Coverage-level thresholding

### The wrong way to pick a contour level

The naive approach is to pick density values directly: "draw a contour at p = 0.05." This is a bad idea for two reasons:

- Density values are not interpretable. `p = 0.05` means nothing to a reader.
- The same density value means different things for different populations. A tight, compact population has high peak density; a diffuse one has low peak density. A fixed density threshold would draw a tiny ring around the tight population and a huge ring around the diffuse one — even though both rings would be at very different *coverage* levels.

### The right way: invert the question

Instead of asking "what region has density above threshold t?", ask "what threshold t gives me a region containing 95% of the mass?" This is the **highest-density region (HDR)** at coverage 95%.

The HDR has a precise, dataset-independent interpretation: "the smallest region containing 95% of this population's probability mass." It is exactly the region you would describe in words as "where this population is."

### The algorithm

Given the KDE and a target coverage c (say 0.95):

1. Evaluate the density at every sample point: `d_i = kde(x_i, y_i)` for each i.
2. Sort those values in descending order. The largest density values correspond to samples in the densest part of the population — the "core."
3. Compute the cumulative sum, normalized: `cdf_k = (d_1 + d_2 + ... + d_k) / sum(d_i)`.
4. Find the smallest k such that `cdf_k >= c`. The density value at that rank, `d_k`, is your threshold.

In code (this is what `points2kde.py` does):

```python
sample_density = kde(sub.T)
sorted_d = np.sort(sample_density)[::-1]      # descending
cdf = np.cumsum(sorted_d) / sorted_d.sum()
thresholds = [sorted_d[np.searchsorted(cdf, c)] for c in [0.95, 0.75, 0.50]]
```

### Why this works

Each sample point is a draw from the true density. The density value at a sample, weighted across all samples, is a Monte-Carlo estimate of the integral of `p^2` — which is what governs HDR computation. Sorting by density and walking the CDF gives you the threshold below which exactly the bottom `(1-c)` fraction of the mass sits. Equivalently, above the threshold is exactly fraction `c` of the mass.

This is the same idea as choosing a quantile, but in density space rather than in any coordinate axis. It is invariant to the shape of the distribution: elongated, multimodal, skewed — the HDR construction handles all of them correctly.

### What you have after step 2

For each population and each coverage level (50%, 75%, 95%) a single density threshold `t`. The level set `{(x, y) : kde(x, y) = t}` is the HDR boundary for that coverage.

---

## 3. Contour rendering

### Evaluate on a grid

A contour-drawing routine like matplotlib's `contour` needs the density evaluated on a regular grid:

```python
xs = np.linspace(x_lo, x_hi, 400)
ys = np.linspace(y_lo, y_hi, 400)
XX, YY = np.meshgrid(xs, ys)
Z = kde(np.vstack([XX.ravel(), YY.ravel()])).reshape(XX.shape)
```

`Z` is a 400×400 array of density values. The grid resolution is a tradeoff between speed (smaller grid = faster) and contour smoothness (larger grid = smoother curves). 400×400 is comfortable for a figure.

### Draw the level sets

Each threshold becomes one `contour` call:

```python
ax.contour(XX, YY, Z, levels=[t95], colors=[color], linewidths=1.8, linestyles='solid')
ax.contour(XX, YY, Z, levels=[t75], colors=[color], linewidths=1.3, linestyles='dashed')
ax.contour(XX, YY, Z, levels=[t50], colors=[color], linewidths=1.0, linestyles='dotted')
```

matplotlib uses the marching squares algorithm internally to trace the level set — for each grid cell it checks which of the four corners are above the threshold and emits line segments accordingly. The result is a closed (or near-closed) curve that follows the density boundary.

Three rings per population, one per coverage level, drawn in the population's color. Solid for the outermost (95%), dashed for the middle (75%), dotted for the innermost (50%) — line style encodes coverage so the color channel is free for population identity.

### Why outline-only and not filled

The first KDE prototype used semi-transparent filled contours. When two populations overlap, the alpha-blend manufactures a *third color* in the overlap region — exactly the failure mode the Figure 1 critique identifies. Outline-only rendering guarantees that **every colored pixel in the chart belongs to exactly one population.** Overlap is conveyed by the geometric intersection of rings, not by any synthesized color.

---

## Putting the pieces together

For each population p:

1. Extract the (x1, x2) coordinates of every sample with label p.
2. Fit `kde_p = gaussian_kde(points_p.T)`.
3. Compute density at each sample: `d_i = kde_p(points_p.T)`.
4. From sorted descending d, find thresholds at 50%, 75%, 95% coverage.
5. Evaluate `kde_p` on the shared (XX, YY) grid.
6. Draw three contours at those three thresholds, in p's color, with solid/dashed/dotted line styles.

Then overlay centroid markers and a thin sample scatter as a background reference, and you have the figure.

---

## What the chart shows, and what it refuses to show

**Shows:**
- The shape of each population's distribution — elongation, skew, multimodality if present
- The relative compactness of populations (tight rings vs. diffuse)
- Where populations overlap — visible as ring intersections
- The position of each population's centroid

**Refuses to show:**
- A categorical label for an arbitrary point in PC space. The chart does not answer "which population is at (x, y)?" because in overlap zones that question has no single correct answer. The honesty of the encoding is in this refusal.

---

## When to reach for this technique

KDE outline contours are the right choice when:

- The labels are categorical but the support is continuous and overlapping.
- The reader needs to see *each population's footprint* and *where footprints overlap*.
- You want a method that makes no parametric assumption about distribution shape.
- You are willing to give up the "color the whole plane" idiom in exchange for not lying about overlap.

They are *not* the right choice when:
- You need a hard assignment for every point — for that, use Voronoi or a soft-assignment GMM as a computational tool (just not as the figure).
- You have so many populations that the rings become visually unreadable. With 20+ populations, color and line-style cannot both stay legible.
- The populations are so well-separated that overlap is not the story — a simpler scatter with centroid markers may communicate more cleanly.
