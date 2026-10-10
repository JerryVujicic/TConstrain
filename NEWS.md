# TConstrain 0.1.0

* Add `fit_tree_posterior()` with a normalized Gaussian centroid likelihood,
  normalized prior on feasible trees, exact posterior edge marginals, conditional
  evidence, entropy and MAP probability. Add `score_tree()` for complete trees.
* Replace the historical `probabilistic = TRUE` distance heuristic with this
  explicit posterior model. This intentionally changes soft-mode results.
* Enforce hard inclusion/exclusion exactly and reject infeasible support.
* Separate absolute-gap temporal supervision (`time_strength`, `time_scale`)
  from pairwise strength. Temporal supervision defaults to zero.
* Collapse duplicate undirected rules, reject conflicts/unknown IDs, preserve
  zero-distance edges, and validate finite coordinates and scalar parameters.
* Support one-dimensional data and singleton clusters. Rooted paths ignore
  optimization weights. Remove direction arrows from undirected rule plots.
* Deprecate `small_eps` and `large_value`; retain the arguments for call compatibility.
* Add structured audits, exhaustive mathematical tests, simulation validation,
  framework documentation and continuous R package checks.
