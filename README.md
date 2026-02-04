⚠️ Note: This package is under development. API is subject to change.

# **TConstrain: Topology-Constrained Trajectory Inference**

**TConstrain** is an R package for semi-supervised single-cell
trajectory inference. It’s a preprocessor-like tool to overcome the
limitations of purely unsupervised methods by allowing researchers to
incorporate biological prior knowledge—such as time points, known
differentiation pathways, or experimental perturbations—as topological
constraints. 

## **Why TConstrain?**

Current state-of-the-art trajectory methods rely on transcriptomic similarity to infer developmental time.
While powerful, this assumption lack flexibility in several biological
scenarios:

**Expression Rebound**: Differentiated cells returning to the original state may be incorrectly mapped as "neighbors" to the start state, creating short-circuits in the trajectory.

**Convergent Differentiation**: Distinct lineages merging into a similar phenotype may be incorrectly linked.

**Disconnected Time-Points**: Large transcriptional jumps between experimental days can cause unsupervised methods to fragment a continuous lineage.

## **How It Works**

TConstrain modifies the Minimum Spanning Tree (MST) construction step
used in standard algorithms. Instead of relying solely on Euclidean
distance, it accepts a user-defined Constraint Table:

**Must-Link**: Forces specific clusters to connect (bridging gaps or enforcing time direction).

**Cannot-Link**: Forbids connections between transcriptionally similar but biologically distinct states (preventing false loops).

**Temporal Penalty**: (Optional) Softly penalizes edges that flow backwards against experimental time labels.

## Installation

You can install the developing version of TConstrain like so:

``` r
# install.packages("devtools")
devtools::install_github("JerryVujicic/TConstrain")
```

## **Potential Future Roadmap**

**Integration with RNA Velocity**: Future versions will support automated constraint suggestion based on scVelo vector fields.

**Probabilistic Constraints**: Implementation of Bayesian priors for edge weights, allowing data-driven overriding of weak biological assumptions.

**Interactive Visualization**: Tools to interactively define and inspect topological anchors on reduced dimension plots.
