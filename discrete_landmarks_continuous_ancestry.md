# Discrete Landmarks, Continuous Ancestry: A Visualization Framework

**Context:** companion to `voronoi_to_kde_decision.md` and to the gamut-clipping analysis in `pca-chart-failure-analysis/pclai_color_ramp_analysis.md`. Those documents make the *negative* case — why the original gradient encoding in PCLAI Figure 1 fails both at the level of cognitive task (asking the viewer to decode a continuous color blend back into ancestry) and at the level of implementation (the CIELAB → sRGB pipeline does not deliver the perceptual uniformity it claims). This document makes the *affirmative* case: what the right encoding should do, and why an encoding built on discrete landmarks with visibly interspersed data is in service of the paper's own thesis rather than at war with it.

The argument is anchored in a claim about the viewer, not about the data: **the human visual system is already capable of performing the inference the original gradient tries to perform for it.** A figure that respects this capability tells a stronger story than a figure that pre-empts it — and, as the implementation analysis shows, a stronger story than a figure whose mathematical foundation has silently collapsed.

---

## 1. The cognitive division of labor

The human visual system is extraordinarily good at certain tasks:

- **Gestalt grouping.** Discrete clusters in 2D are recognized as named entities immediately, without conscious effort. The visual system organizes the field into figures and grounds before any deliberate reading happens.
- **Inter-landmark interpolation.** A point that lies between two named clusters is perceived *as lying between them*. The viewer reads spatial intermediacy directly, the same way they read "halfway between two mountains" on a map.
- **Boundary-crossing detection.** When a point sits across or outside a named region, the viewer notices. The cognitive event is automatic: *this one is not where the others are.*

The human visual system is extraordinarily bad at certain other tasks:

- Decoding precise sub-JND color differences across a continuous field.
- Reverse-engineering "what proportion of red, green, blue went into this muddy beige?"
- Translating a position in a continuous color gradient back into a quantitative claim about the underlying weights.

The gamut-clipping analysis in `pca-chart-failure-analysis/pclai_color_ramp_analysis.md` shows that the original Figure 1 background asks the viewer to do the things they are bad at — and, worse, asks them to do so against a color field that does not deliver what it advertises. The stated design rationale is CIELAB perceptual uniformity: equal distances in (whitened) PCA space should produce equal perceived color differences, so that color similarity tracks genetic similarity. That guarantee holds only *within* the sRGB gamut. The implementation maps the reference data to `a*, b* ∈ [−120, +120]` at a fixed `L* = 80`, where the sRGB gamut occupies only **16.4%** of the declared bounding box. As a consequence, **100% of the 3,122 reference haplotypes fall outside the gamut** and are clipped to its boundary — a closed 1D curve. The intended 2D perceptually-uniform color surface collapses to a 1D contour, the clipping is unequal across continental groups (European haplotypes displaced 69 Lab units on average, South Asian 21), and the final figure is then alpha-blended at 0.38 onto white, desaturating the already-distorted colors a second time.

This document makes the complementary affirmative point: **the same data, encoded with discrete landmarks and a visible scatter, asks the viewer to do the things they are good at — and does so with an encoding whose visual claims it can actually keep.**

The encoding should get out of the way of capacities the viewer already possesses, not substitute machine-computed interpolation for perceptual interpolation. The perceptual one is, in this domain, more accurate — because it is conditioned on context (relative cluster sizes, shape, density, scatter, neighborhood) that no inverse-distance kernel sees, and because it does not depend on a color-space round-trip that silently fails.

---

## 2. Why this matters for the PCLAI thesis specifically

The paper's central claim is that ancestry along the genome lives at *continuous coordinates*, not at discrete population labels. The whole methodological contribution is the move away from forced categorical assignment.

A continuous gradient background visually argues the opposite thesis. By coloring every point in the (PC1, PC2) plane with a deterministic ancestry blend, it asserts that every coordinate has a definite, computable categorical interpretation — that there is no genuine "between." That is exactly the position the paper is rejecting in its text.

The original encoding is at war with the paper's own argument — and, separately, with itself. The CIELAB choice was made to honor a continuous-distance reading of ancestry: equal color differences should mean equal genetic distances. The implementation does not preserve that property for a single point actually displayed in the figure. The gradient simultaneously over-claims (treating intermediate coordinates as if they had categorical answers) and under-delivers (the color that encodes those answers is a clipped 1D projection, not the intended 2D perceptual surface). Both failures point the same direction: the figure is doing work the encoding cannot honestly support.

A discrete-landmark encoding flips this. It shows:

- The reference populations as **fixed historical landmarks** — geographically and anthropologically meaningful entities with bounded footprints in PC space.
- The data — including admixed individuals and, in the methods context, haplotype-level segment coordinates — as **points that may live inside one footprint, between several, or outside them all.**

The viewer sees, directly and immediately, the phenomenon the paper is about: *coordinates that do not reduce to categories.* The chart becomes an illustration of the thesis rather than a contradiction of it.

The intermediate points are not noise to be smoothed away by a background gradient. They are the finding.

---

## 3. The agency principle

A useful way to state the principle behind the encoding:

> The figure should provide the **raw material** for inference, not the **conclusion** of inference.

The viewer is performing a scientific act when they look at the chart. They are forming a judgment about where a point sits relative to named reference populations and what that location means. An encoding that pre-computes that judgment — whether by Voronoi cells or by a continuous color field — strips the viewer of the act they came to the figure to perform.

The discrete landmarks supply what the viewer cannot derive: the historical and anthropological identity of each named reference population, computed once from labeled samples, drawn cleanly on the plane. The visible scatter supplies what the viewer is fully equipped to interpret: the spatial relationships between data and landmarks. The encoding stops there. Everything else — does this point belong, is it intermediate, is it ambiguous — is the viewer's call.

This is the same principle that makes Tufte's "show the data" injunction work: maximize the information density of the things the viewer can actually use, minimize the pre-digested judgment baked into the encoding.

---

## 4. KDE outline contours as the realization

The algorithm and rendering details are in `kde_tutorial.md`. What matters for the framing here is that the technique has exactly the properties the principle requires — and that it sidesteps the specific failure mode the gamut analysis exposes:

- **Each landmark is a distinct visual entity.** Three nested rings at 50% / 75% / 95% coverage define the footprint of each reference population. The rings have a precise verbal interpretation (HDR coverage), so the named landmark is not a vague halo — it is a quantitatively defined region.
- **Every colored pixel belongs to exactly one population.** No synthesized colors appear in regions of shared support; overlap is shown by ring intersection. The encoding refuses to manufacture an answer where the data has multiple answers.
- **The scatter is preserved.** Light per-population sample points sit beneath the rings, so the viewer can always see the underlying data — not just the smoothed footprint, but the actual individuals contributing to it.
- **Population shape is preserved.** Elongation, skew, multimodality within a single population are visible in the ring geometry. The encoding does not impose an isotropic or Gaussian shape on what may be neither.
- **No load is placed on perceptual color-distance.** Color in this encoding is used only as a categorical landmark label, drawn from a small, well-separated, in-gamut palette. The encoding does not claim that color similarity tracks genetic distance, so it cannot fail to deliver that claim. The quantitative content lives in geometry (ring position, ring overlap, point location relative to rings), where the perceptual system is reliable, rather than in sub-JND color differences inside a continuous field, where it is not — and where, in the original implementation, the underlying field is in any case clipped to the gamut boundary.

The combination produces the cognitive setup the principle calls for: discrete landmarks the viewer can identify at a glance, and visible data the viewer can locate relative to them.

---

## 5. A federation of visualizations

There is no one true visualization. A figure is a compression of the data for a specific audience and a specific question, and different questions are served by different compressions. The right stance is **a federation of visualizations** — multiple complementary views, each making explicit what the others leave implicit.

Within this federation, the KDE outline encoding is the primary view: it serves the paper's thesis directly, supports the viewer's agency, and stays within the perceptual envelope the four-layer critique demands.

A Voronoi tessellation of population centroids has a legitimate place in the federation — not as a competing primary view, but as a **historical / methodological companion**:

- **Historical reference.** The Voronoi partition expresses the legacy assumption that ancestry is a categorical assignment problem. Including it makes that assumption visible rather than invisible, and lets the reader see how the field traditionally framed the question the paper is now reframing.
- **Algorithmic contrast.** For readers comparing the PCLAI approach to nearest-centroid or hard-assignment baselines, the Voronoi panel is the natural visual proxy for that baseline. The reader can see, directly, what is lost when categorical assignment is forced onto a continuous structure.
- **Decision boundary intuition.** Even as a tool the paper rejects for primary visualization, Voronoi remains a useful pedagogical device: every reader knows what a nearest-centroid partition is. Showing it explicitly defangs it.

A composite figure pairing the two — Voronoi on one panel, KDE on another, same data, same coordinate system — would let the reader form the comparison the paper's text is making. This is stronger than either chart alone: the KDE shows what the paper *is* doing; the Voronoi shows what the paper is *not* doing. The federation makes the move visible.

Other federation members worth considering:

- A small-multiples panel of (PC1, PC3), (PC2, PC3), etc., to recover structure that any single 2D projection loses.
- A purely categorical scatter (sample points colored by superpopulation, no background field at all) as a minimal reference — useful as a "trust no synthesis" view.

The principle is that each view answers a specific question well, and that the figure earns its complexity by serving multiple questions rather than by amassing visual layers in service of one.

---

## 6. What the viewer is meant to do

A final way to test an encoding is to ask: what cognitive act does it invite?

| Encoding | What it invites the viewer to do |
|---|---|
| Continuous CIELAB gradient (as implemented) | Decode a sub-JND color blend into a categorical ancestry guess — from a field that is, in fact, a 1D gamut-boundary projection of the intended 2D surface, unequally distorted across continental groups. (Bad at this; impossible in this implementation.) |
| Voronoi tessellation | Read a hard-boundary categorical label off the partition. (Wrong question.) |
| KDE outline contours | Locate a point relative to named footprints and judge intermediacy. (Good at this; right question; no perceptual-uniformity claim to break.) |

The KDE encoding invites exactly the cognitive act the paper is asking the reader to perform: see the discrete landmarks, see the data, see the continuous coordinate, decide for yourself what intermediacy means.

That is the figure the paper deserves.

---

## Summary

The argument for the KDE outline encoding is not merely defensive (it avoids dishonesty) or technical (it has the right perceptual properties). The strongest case is constructive: the encoding **serves the paper's own thesis** by giving the viewer the raw material for a continuous-coordinate reading of ancestry, rather than pre-empting that reading with a categorical color field.

It also avoids a specific, quantifiable failure that the current implementation has not escaped. The CIELAB gradient was meant to encode genetic distance as perceptual distance; the gamut analysis shows that 100% of reference haplotypes are clipped to the sRGB boundary, collapsing the intended 2D surface to a 1D contour and introducing systematic, group-dependent distortions. The KDE encoding makes no such promise to break: it carries quantitative content in geometry, uses color only as a categorical landmark label drawn from an in-gamut palette, and leaves the continuous-coordinate reading to the viewer — who is, in this domain, the more reliable instrument.

The encoding works because the human visual system already does the interpretive work that matters. The right job of the figure is to lay down clean, named landmarks; to show the actual data points; and then to step aside. The Voronoi tessellation belongs in the federation as a historical companion that makes the rejected baseline explicit, but the primary view should be the one that lets continuous ancestry look like continuous ancestry — without quietly subcontracting that job to a color pipeline that cannot perform it.
