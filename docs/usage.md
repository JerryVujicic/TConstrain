# TConstrain usage guide

[Back to the README](../README.md)

This guide covers configuration, interpretation and migration. For the model
equations and derivation, see the
[probabilistic framework](../vignettes/probabilistic_framework.Rmd).

## Inputs

`fit_tree_posterior()` accepts a finite numeric matrix or data frame of cluster
centroids, with one cluster per row and unique cluster IDs as row names.
`get_centroids(data, labels)` calculates these from cell coordinates and matching
cluster labels. `infer_trajectory()` combines centroid calculation, tree inference
and rooted cluster paths.

Constraints are data frames with `from`, `to` and `type` columns. Both endpoints
must be known cluster IDs, and `type` is either `must_link` or `cannot_link`.
These rules concern direct, undirected adjacency. A cannot-link still permits
an indirect path between its endpoints.

## Soft constraints and parameter choice

```r
library(TConstrain)
centroids <- rbind(A = c(0, 0), B = c(1, 1),
                  C = c(0.2, 0), D = c(2, 1))
soft <- data.frame(from = "A", to = "C", type = "cannot_link")
fit <- fit_tree_posterior(centroids, constraints = soft,
                          likelihood_scale = 1, prior_strength = 5)
```

The Gaussian likelihood models centroid coordinates as independent isotropic
increments along a proposed tree, conditional on one fixed centroid. The
anchored centroid is the first row. Changing the reporting root does not supply
directional evidence. The model introduces no latent ancestral nodes and does
not account for centroid estimation uncertainty.

`prior_strength` is a prior weight ratio for a one-edge exchange against a
neutral edge. It is not the marginal probability that a link exists. A soft
must-link favors an edge; a soft cannot-link discourages it. Neither guarantees
an outcome. The prior is normalized over all spanning trees satisfying the hard
constraints, and the likelihood updates it to a normalized posterior.

`likelihood_scale` defaults to **1 coordinate unit**, not a learned or calibrated
value. Rescale it with the coordinates. Set it from external information or
assess a grid on held-out biological evidence. The package does not infer
hyperparameters or integrate over their uncertainty. Choosing settings using
the same fitted evidence is an empirical tuning step, not full Bayes.

## Exact hard constraints

Hard rules can be combined with soft preferences in a single fit:

```r
hard <- data.frame(from = "A", to = "D", type = "cannot_link")
fit <- fit_tree_posterior(centroids, constraints = soft,
                          hard_constraints = hard, likelihood_scale = 1,
                          prior_strength = 5)
```

Required edges are contracted before normalization/optimization; excluded edges
are removed. Cycles, conflicting rules and disconnected support raise errors.
No finite "expensive bridge" overrides a hard exclusion. A disconnected system
should be analyzed separately; this release does not fit distributions on forests.
Hard rules dominate soft preferences on the same pair.

Duplicate undirected rules are collapsed with a warning; unknown IDs fail.
Zero-distance edges and single-cluster/one-dimensional inputs are supported.

For a geometric constrained MST without a posterior, use
`get_constrained_mst(centroids, constraints = hard)`. With its default
`probabilistic = FALSE`, the supplied constraints are hard and the objective is
Euclidean distance plus any independently specified temporal gap penalty.

## Rooted cluster paths

For cell-level inputs, the wrapper returns ordered cluster paths:

```r
# A small demonstration; replace these repeated coordinates with your cells.
coords <- centroids[rep(seq_len(nrow(centroids)), each = 2), ]
labels <- rep(rownames(centroids), each = 2)
trajectory <- infer_trajectory(coords, labels, start_cluster = "A",
    probabilistic = TRUE, constraints = soft, hard_constraints = hard,
    likelihood_scale = 1, prior_strength = 5)
trajectory$lineages   # Ordered cluster names, not cell pseudotimes
trajectory$posterior # Full posterior result
```

Specify `start_cluster`, or supply time labels with a unique earliest known
cluster. An explicit root overrides time-based selection. The root sets the
reference for reading paths; it does not turn undirected edges into inferred
biological transitions.

## Time labels

Time supervision is **off by default** (`time_strength = 0`) and is separate
from pairwise strength. Positive strength favors smaller absolute gaps in the
units set by `time_scale`. `NA` means unknown. This is an undirected prior, not a
chronological ordering guarantee. Time labels can still choose a unique earliest
root even when temporal supervision is off.

```r
# Illustrative times and strengths, not fitted biological parameters.
times <- c(A = 0, B = 1, C = 2, D = 3)
timed_fit <- fit_tree_posterior(centroids, constraints = soft,
    time_labels = times, likelihood_scale = 1, prior_strength = 5,
    time_strength = 1, time_scale = 1)
```

## Posterior results and audit

`fit$edges` exposes original distances, pair/time log prior contributions,
likelihoods, prior/posterior marginals, hard support and MAP membership.
Edges are dependent because they must form a tree: their marginal probabilities
sum to `number of clusters - 1`, not 1. Their product is not a tree probability.

```r
fit$edges[, c("from", "to", "prior_probability", "posterior_probability", "in_map")]
fit$map_probability
score_tree(fit, fit$map)
```

`fit$map_probability` is the posterior mass of the reported complete tree.
`score_tree()` calculates the likelihood, normalized prior and posterior for
a candidate tree. `fit$entropy` describes uncertainty across trees in nats.

`fit$log_evidence` is a conditional coordinate log density. On the density
scale, evidence can exceed one. Compare it only for the same coordinates, units,
dimensions, anchor and observation model.

Both inference modes return structured audits. Edge-list returns from
`get_constrained_mst()` preserve audit, posterior and vertex IDs as R attributes.
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

See the [release notes](../NEWS.md) for the complete change summary.

## Validation and scientific scope

Run `testthat::test_local()` for exhaustive small-graph comparisons and regression
tests. `Rscript validation/run_validation.R` reproduces the rebound analysis,
prior sensitivity and a Gaussian-model simulation study. The
[validation report](../validation/REPORT.md) records the results and limitations.

The output is cluster-level topology, not cell pseudotime, directed transitions,
ancestry, an RNA-count model or fate probabilities. A tree cannot describe GC
cycling or multi-parent convergence. Use the results as inspectable topology
hypotheses and evaluate them with independent biological evidence.

The [framework vignette](../vignettes/probabilistic_framework.Rmd) gives the
likelihood and normalization equations, numerical safeguards, complexity and
method references. Its Gaussian coordinate model is a specific observation
model, distinct from the categorical belief-network model in the cited work.
