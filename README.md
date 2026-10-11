# TConstrain

TConstrain infers **cluster-level trajectory trees** from reduced single-cell
coordinates. It supports exact must-link/cannot-link constraints and
probabilistic inference with a Gaussian likelihood, normalized tree priors,
and posterior edge probabilities.

The output is an undirected cluster topology; cell-level pseudotime and fate
probabilities require additional methods.

## Installation

```r
# install.packages("remotes")
remotes::install_github("JerryVujicic/TConstrain")
```

## Quick start

```r
library(TConstrain)

# Illustrative cluster coordinates and a soft adjacency constraint.
centroids <- rbind(A = c(0, 0), B = c(1, 1),
                  C = c(0.2, 0), D = c(2, 1))
soft <- data.frame(from = "A", to = "C", type = "cannot_link")

fit <- fit_tree_posterior(centroids, constraints = soft,
                          likelihood_scale = 1, prior_strength = 5)
fit$map  # Maximum a posteriori tree
fit$edges[, c("from", "to", "posterior_probability")]
```

`likelihood_scale` is a fixed standard deviation in coordinate units, not an
estimated parameter. See the usage guide for choosing settings and interpreting
probabilities.

## Documentation

- [Usage guide](docs/usage.md): inputs, constraints, time labels, outputs and migration.
- [Probabilistic framework](vignettes/probabilistic_framework.Rmd): likelihood, prior, posterior and assumptions.
- [Rebound example](vignettes/rebound_demo.Rmd): an illustrated trajectory analysis.
- [Validation report](validation/REPORT.md): tests, simulations and remaining limitations.
- [Release notes](NEWS.md): changes in version 0.1.0.
