# TConstrain

TConstrain infers **cluster-level undirected trajectory trees** from reduced
single-cell coordinates. Version 0.1.0 adds a normalized tree prior, an explicit
Gaussian likelihood, exact posterior edge probabilities and a MAP tree. It also
enforces hard must-link and cannot-link rules literally.

The package does not calculate cell pseudotime, directed transitions, ancestry,
RNA-count likelihoods or fate probabilities. A tree cannot describe GC cycling
or multi-parent convergence. Use the results as inspectable topology hypotheses.

## Installation

```r
# install.packages("remotes")
remotes::install_github("JerryVujicic/TConstrain")
```

## Probabilistic inference

```r
library(TConstrain)
centroids <- rbind(Naive = c(0, 0), Activated = c(1, 1),
                  Memory = c(0.2, 0), Plasma = c(2, 1))
# Illustrative hypotheses, not established B-cell lineage relationships.
soft <- data.frame(from = "Naive", to = "Memory", type = "cannot_link")
fit <- fit_tree_posterior(centroids, constraints = soft,
                          likelihood_scale = 1, prior_strength = 5)
fit$edges[, c("from", "to", "prior_probability", "posterior_probability", "in_map")]
fit$map_probability
score_tree(fit, fit$map)
```

For coordinates `x` in `d` dimensions and a fixed increment standard deviation
`sigma = likelihood_scale`, the likelihood conditional on one fixed centroid is

```
p(X[-r] | X[r], T, sigma) = product over (i,j) in T of
    (2*pi*sigma^2)^(-d/2) * exp(-||x[i]-x[j]||^2/(2*sigma^2)).
```

This is a normalized Gaussian random-walk density on the observed centroid
nodes. It assumes independent isotropic increments along a proposed tree and
conditions on coordinates as observed, without centroid estimation uncertainty.
The anchored centroid is the first row; changing the reporting root does not
supply directional evidence. No ancestral latent nodes are introduced.

Let `a[i,j]` equal `+log(prior_strength)` for a soft must-link,
`-log(prior_strength)` for a soft cannot-link, and zero otherwise. Subtract
`time_strength * abs(time[i]-time[j]) / time_scale` for two known times.
The tree prior is `exp(sum(a[T])) / Z_prior`, on **all spanning trees satisfying
hard constraints**. The posterior is proportional to this prior times the
likelihood, with its own normalizer. The matrix-tree theorem computes both
normalizers and exact edge marginals without enumerating trees.

`prior_strength` is a prior weight ratio for a one-edge exchange against a
neutral edge. It is not the marginal probability that a link exists. Edges
are dependent because they must form a tree; their marginal probabilities sum
to `number of clusters - 1`, not 1.

`likelihood_scale` defaults to **1 coordinate unit**, not a learned or calibrated
value. Rescale it with the coordinates. Set it from external information or
assess a grid on held-out biological evidence. The package does not infer
hyperparameters or integrate over their uncertainty. Choosing settings using
the same fitted evidence is an empirical tuning step, not full Bayes.

## Hard rules and rooted paths

```r
hard <- data.frame(from = "Naive", to = "Plasma", type = "cannot_link")
fit <- fit_tree_posterior(centroids, constraints = soft, hard_constraints = hard)

# coords: cells x dimensions; labels: one cluster ID per cell
# trajectory <- infer_trajectory(coords, labels, start_cluster = "Naive",
#     probabilistic = TRUE, constraints = soft, hard_constraints = hard,
#     likelihood_scale = 1, prior_strength = 5)
# trajectory$lineages   # Ordered cluster names, not cell pseudotimes
# trajectory$posterior # Full posterior result
```

Required edges are contracted before normalization/optimization; excluded edges
are removed. Cycles, conflicting rules and disconnected support raise errors.
No finite "expensive bridge" overrides a hard exclusion. A disconnected system
should be analyzed separately; this release does not fit distributions on forests.
Duplicate undirected rules are collapsed with a warning; unknown IDs fail.
Zero-distance edges and single-cluster/one-dimensional inputs are supported.

## Time and audit

Time supervision is **off by default** (`time_strength = 0`) and is separate
from pairwise strength. Positive strength favors smaller absolute gaps in the
units set by `time_scale`. NA is unknown. This is an undirected prior, not a
chronological ordering guarantee. Time labels can still choose a unique earliest
root even when temporal supervision is off. An explicit root overrides it.

`fit$edges` exposes original distances, pair/time log prior contributions,
likelihoods, prior/posterior marginals, hard support and MAP membership.
`fit$log_evidence` is a conditional coordinate density, and can exceed one on
the density scale. Compare it only for the same coordinates, units, dimensions,
anchor and observation model. `fit$entropy` describes tree uncertainty in nats.

Hard MST mode (`probabilistic = FALSE`) remains available. Its cost is distance
plus the independent temporal gap penalty. Both modes return structured audits;
edge-list returns preserve audit, posterior and vertex IDs as R attributes.
MAP costs can be negative; they must never be read as developmental duration.

## Migration from the development version

- `probabilistic = TRUE` now uses the Gaussian posterior MAP model. It no longer
  divides or multiplies geometric distances. Old strength values need reassessment.
- Time labels alone no longer apply a strength-coupled same/different multiplier.
  Set `time_strength` and `time_scale` explicitly when temporal supervision is wanted.
- `small_eps` and `large_value` are deprecated and ignored with a warning when
  explicitly supplied. Hard constraints no longer depend on data scale.
- Invalid or contradictory rules error instead of being silently dropped.
- Zero distances remain edges; paths use unweighted graph connectivity.

## Validation and limits

Run `testthat::test_local()` for exhaustive small-graph comparisons and regression
tests. `Rscript validation/run_validation.R` reproduces the rebound analysis,
prior sensitivity and a Gaussian-model simulation study. See
[the validation report](validation/REPORT.md) and
[the framework vignette](vignettes/probabilistic_framework.Rmd).
GitHub Actions runs R package checks on pushes and pull requests.

The dense implementation uses O(K^2) memory and O(K^3) linear algebra for K
clusters. Underflow or ill-conditioned normalizers cause informative errors,
not silently truncated probabilities. For uncertain or poorly specified
geometry, inspect sensitivity across embeddings, clusters, donors and model
parameters. Posterior certainty under this model does not establish biological
correctness. Shared BCR clonotype is not a direct adjacency observation; reserve
independent lineage/time/functional measurements for evaluation.

The factorized-tree normalization follows the matrix-tree framework described
by [Meila and Jaakkola (UAI 2000)](https://arxiv.org/abs/1301.3875).
The Gaussian coordinate model here is a specific, deliberately limited
observation model, not their categorical belief-network model.
