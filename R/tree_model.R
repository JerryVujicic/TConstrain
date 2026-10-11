# Internal validation and tree algorithms shared by both inference modes.
.tc_scalar <- function(x, name, lower = 0, strict = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      if (strict) x <= lower else x < lower) {
    stop(sprintf("'%s' must be a finite numeric scalar %s %s.",
                 name, if (strict) ">" else ">=", lower), call. = FALSE)
  }
}

.tc_coordinates <- function(x, named = TRUE) {
  if (!is.matrix(x) && !is.data.frame(x)) stop("Coordinates must be a matrix or data.frame.")
  x <- as.matrix(x)
  if (!is.numeric(x) || !nrow(x) || !ncol(x) || any(!is.finite(x))) {
    stop("Coordinates must be nonempty, numeric and finite.")
  }
  if (named) .tc_ids(rownames(x), nrow(x))
  x
}

.tc_ids <- function(ids, n) {
  if (is.null(ids) || length(ids) != n || anyNA(ids) ||
      any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("Cluster IDs must be unique, nonmissing and nonempty.")
  }
}

.tc_rules <- function(x, ids) {
  if (is.null(x)) return(data.frame(i = integer(), j = integer(), type = character()))
  if (!is.data.frame(x) || !all(c("from", "to", "type") %in% names(x))) {
    stop("Constraints must be a data.frame with 'from', 'to' and 'type'.")
  }
  a <- match(as.character(x$from), ids)
  b <- match(as.character(x$to), ids)
  type <- as.character(x$type)
  if (anyNA(a) || anyNA(b)) stop("Constraints contain missing or unknown cluster IDs.")
  if (anyNA(type) || any(!type %in% c("must_link", "cannot_link"))) stop("Unknown constraint type.")
  if (any(a == b)) stop("Self constraints are not allowed.")
  out <- data.frame(i = pmin(a, b), j = pmax(a, b), type = type)
  key <- paste(out$i, out$j, sep = ":")
  if (any(vapply(split(type, key), function(z) length(unique(z)) > 1L, logical(1)))) {
    stop("Conflicting constraints on the same undirected pair.")
  }
  if (anyDuplicated(key)) warning("Duplicate undirected constraints collapsed to one rule.", call. = FALSE)
  out[!duplicated(key), , drop = FALSE]
}

.tc_times <- function(time_labels, ids) {
  if (is.null(time_labels)) return(rep(NA_real_, length(ids)))
  if (!is.numeric(time_labels) || is.null(names(time_labels))) stop("'time_labels' must be a named numeric vector.")
  .tc_ids(names(time_labels), length(time_labels))
  if (any(is.infinite(time_labels)) || any(is.nan(time_labels))) stop("Times must be finite or NA.")
  if (any(!names(time_labels) %in% ids)) stop("'time_labels' contains unknown cluster IDs.")
  unname(time_labels[ids])
}

.tc_edges <- function(x, time_labels, time_scale) {
  n <- nrow(x)
  ij <- which(upper.tri(matrix(FALSE, n, n)), arr.ind = TRUE)
  d <- as.matrix(stats::dist(x))
  if (any(!is.finite(d))) stop("Pairwise distances overflowed; rescale the coordinates.")
  times <- .tc_times(time_labels, rownames(x))
  gap <- abs(times[ij[, 1]] / time_scale - times[ij[, 2]] / time_scale)
  gap[is.na(gap)] <- 0
  if (any(!is.finite(gap))) stop("Time gaps overflowed; rescale the times.")
  data.frame(from = rownames(x)[ij[, 1]], to = rownames(x)[ij[, 2]],
             i = ij[, 1], j = ij[, 2], distance = d[ij], time_gap = gap)
}

.tc_rule_indices <- function(rules, edges) {
  match(paste(rules$i, rules$j, sep = ":"), paste(edges$i, edges$j, sep = ":"))
}

.tc_support <- function(n, edges, rules) {
  required <- forbidden <- rep(FALSE, nrow(edges))
  idx <- .tc_rule_indices(rules, edges)
  required[idx[rules$type == "must_link"]] <- TRUE
  forbidden[idx[rules$type == "cannot_link"]] <- TRUE
  component <- seq_len(n)
  for (k in which(required)) {
    a <- component[edges$i[k]]
    b <- component[edges$j[k]]
    if (a == b) stop("Hard must-link constraints contain a cycle; no spanning tree is feasible.")
    component[component == b] <- a
  }
  component <- match(component, unique(component))
  eligible <- !forbidden & !required & component[edges$i] != component[edges$j]
  # Feasibility check also catches a two-node hard cannot-link.
  connected <- component
  for (k in which(eligible)) {
    a <- connected[edges$i[k]]
    b <- connected[edges$j[k]]
    connected[connected == b] <- a
  }
  if (length(unique(connected)) != 1L) {
    stop("Hard constraints disconnect the graph; no spanning tree exists. Analyze components separately.")
  }
  list(required = required, forbidden = forbidden, eligible = eligible,
       component = component, n_components = max(component))
}

.tc_map <- function(n, edges, support, cost) {
  if (any(!is.finite(cost))) stop("Non-finite edge scores; rescale model inputs.")
  selected <- which(support$required)
  component <- support$component
  candidates <- which(support$eligible)
  candidates <- candidates[order(cost[candidates], edges$i[candidates], edges$j[candidates])]
  for (k in candidates) {
    a <- component[edges$i[k]]
    b <- component[edges$j[k]]
    if (a != b) {
      selected <- c(selected, k)
      component[component == b] <- a
    }
    if (length(selected) == n - 1L) break
  }
  sort(selected)
}

.tc_graph <- function(ids, edges, selected, cost) {
  tab <- edges[selected, , drop = FALSE]
  tab$weight <- cost[selected]
  igraph::graph_from_data_frame(tab, directed = FALSE, vertices = data.frame(name = ids))
}

.tc_matrix <- function(ids, edges, values, diagonal = 0) {
  out <- matrix(diagonal, length(ids), length(ids), dimnames = list(ids, ids))
  out[cbind(edges$i, edges$j)] <- values
  out[cbind(edges$j, edges$i)] <- values
  out
}

# Matrix-tree normalizer after contraction of the mandatory forest. Parallel
# edges are retained in the incidence matrix, so each original edge has its
# own marginal. No probabilities are approximated by thresholding weights.
.tc_tree_distribution <- function(edges, support, log_weights) {
  if (any(!is.finite(log_weights))) stop("Non-finite log weights; rescale model inputs.")
  marginal <- as.numeric(support$required)
  forced <- sum(log_weights[support$required])
  m <- support$n_components
  if (m == 1L) return(list(log_z = forced, marginal = marginal, rcond = 1, entropy = 0))
  use <- which(support$eligible)
  shift <- max(log_weights[use])
  weights <- exp(log_weights[use] - shift)
  if (any(weights == 0)) {
    stop("Tree weights underflowed; increase likelihood_scale or reduce prior/time strength.")
  }
  a <- support$component[edges$i[use]]
  b <- support$component[edges$j[use]]
  # Accumulate parallel conductances without an edges-by-vertices matrix.
  laplacian <- matrix(0, m, m)
  for (k in seq_along(use)) {
    u <- a[k]
    v <- b[k]
    laplacian[u, u] <- laplacian[u, u] + weights[k]
    laplacian[v, v] <- laplacian[v, v] + weights[k]
    laplacian[u, v] <- laplacian[u, v] - weights[k]
    laplacian[v, u] <- laplacian[v, u] - weights[k]
  }
  reduced <- laplacian[-m, -m, drop = FALSE]
  conditioning <- rcond(reduced)
  if (!is.finite(conditioning) || conditioning < 1e-12) {
    stop("Tree normalizer is numerically ill-conditioned; increase likelihood_scale or reduce prior/time strength.")
  }
  ch <- chol(reduced)
  log_z <- forced + (m - 1) * shift + 2 * sum(log(diag(ch)))
  inverse <- matrix(0, m, m)
  inverse[-m, -m] <- chol2inv(ch)
  probability <- weights * (diag(inverse)[a] + diag(inverse)[b] -
                              2 * inverse[cbind(a, b)])
  if (any(probability < -1e-7 | probability > 1 + 1e-7) ||
      abs(sum(probability) - (m - 1)) > 1e-6 * m) {
    stop("Numerical edge-marginal check failed; rescale model inputs.")
  }
  marginal[use] <- pmax(0, pmin(1, probability))
  list(log_z = log_z, marginal = marginal, rcond = conditioning,
       entropy = max(0, log_z - sum(marginal * log_weights)))
}

#' Fit a normalized posterior distribution over cluster trees
#'
#' Conditional on one fixed centroid, coordinates follow independent isotropic
#' Gaussian increments along the edges of a tree. All trees have the same number
#' of increments. The likelihood includes their Gaussian density constants.
#' A normalized edge-factorized tree prior encodes soft rules and time gaps.
#' Hard rules restrict the support exactly. See the probabilistic framework
#' vignette for equations, assumptions and interpretation.
#'
#' @param centroids Numeric matrix or data.frame with unique cluster row names.
#' @param constraints Optional soft rules with columns from, to, type.
#'   Identical undirected duplicates are collapsed; conflicting rules are errors.
#' @param hard_constraints Optional rules with the same columns, enforced exactly.
#'   A cyclic must-link set or disconnected support is an error.
#' @param time_labels Named numeric sampling times; NA denotes unknown time.
#' @param prior_strength Finite number >= 1. A soft must-link contributes
#'   log(prior_strength), and a cannot-link its negative, to an edge log weight.
#'   This is a tree weight ratio for a one-edge exchange, not an edge probability.
#' @param likelihood_scale Finite positive Gaussian increment standard deviation
#'   in coordinate units. Default 1 is fixed, not estimated or calibrated.
#' @param time_strength Finite nonnegative coefficient of absolute time gap in
#'   the log prior. Default 0 disables time supervision independently of rules.
#' @param time_scale Finite positive divisor converting sampling-time units.
#' @return A list of class tconstrain_posterior containing map (an igraph tree),
#'   edges (likelihood, prior and posterior audit), log_prior_normalizer,
#'   log_posterior_normalizer, log_evidence, map_log_probability,
#'   map_probability, entropy, parameters and numerical diagnostics.
#'   Probabilities are conditional on the supplied centroids and fixed model
#'   parameters; they are not cell fate or directed transition probabilities.
#' @examples
#' x <- rbind(A = c(0, 0), B = c(1, 0), C = c(0, 1))
#' fit <- fit_tree_posterior(x, likelihood_scale = 1)
#' fit$edges[, c("from", "to", "posterior_probability")]
#' score_tree(fit, fit$map)
#' @export
fit_tree_posterior <- function(centroids, constraints = NULL, hard_constraints = NULL,
                               time_labels = NULL, prior_strength = 10,
                               likelihood_scale = 1, time_strength = 0, time_scale = 1) {
  x <- .tc_coordinates(centroids)
  ids <- rownames(x)
  .tc_scalar(prior_strength, "prior_strength", 1, FALSE)
  .tc_scalar(likelihood_scale, "likelihood_scale")
  .tc_scalar(time_strength, "time_strength", 0, FALSE)
  .tc_scalar(time_scale, "time_scale")
  edges <- .tc_edges(x, time_labels, time_scale)
  rules <- .tc_rules(constraints, ids)
  hard <- .tc_rules(hard_constraints, ids)
  support <- .tc_support(nrow(x), edges, hard)
  edges$pair_log_prior <- rep(0, nrow(edges))
  idx <- .tc_rule_indices(rules, edges)
  edges$pair_log_prior[idx] <- ifelse(rules$type == "must_link", 1, -1) * log(prior_strength)
  edges$time_log_prior <- -time_strength * edges$time_gap
  edges$log_prior_weight <- edges$pair_log_prior + edges$time_log_prior
  edges$log_likelihood <- -ncol(x) * (log(2 * pi) / 2 + log(likelihood_scale)) -
    (edges$distance / likelihood_scale)^2 / 2
  edges$log_posterior_weight <- edges$log_prior_weight + edges$log_likelihood
  prior <- .tc_tree_distribution(edges, support, edges$log_prior_weight)
  posterior <- .tc_tree_distribution(edges, support, edges$log_posterior_weight)
  selected <- .tc_map(nrow(x), edges, support, -edges$log_posterior_weight)
  edges$required <- support$required
  edges$forbidden <- support$forbidden
  edges$eligible <- support$eligible
  edges$prior_probability <- prior$marginal
  edges$posterior_probability <- posterior$marginal
  edges$in_map <- seq_len(nrow(edges)) %in% selected
  log_p <- sum(edges$log_posterior_weight[selected]) - posterior$log_z
  if (log_p > 1e-7) stop("Numerical MAP probability check failed.")
  log_p <- min(0, log_p)
  ans <- list(map = .tc_graph(ids, edges, selected, -edges$log_posterior_weight),
              edges = edges, cluster_ids = ids,
              log_prior_normalizer = prior$log_z,
              log_posterior_normalizer = posterior$log_z,
              log_evidence = posterior$log_z - prior$log_z,
              map_log_probability = log_p, map_probability = exp(log_p),
              entropy = posterior$entropy,
              parameters = list(likelihood_scale = likelihood_scale,
                                prior_strength = prior_strength,
                                time_strength = time_strength, time_scale = time_scale,
                                dimensions = ncol(x), likelihood_anchor = ids[1]),
              diagnostics = list(prior_rcond = prior$rcond, posterior_rcond = posterior$rcond,
                                 expected_edges = sum(posterior$marginal)))
  class(ans) <- "tconstrain_posterior"
  ans
}

#' Score a complete candidate tree under a fitted model
#'
#' @param fit A result from fit_tree_posterior.
#' @param tree An igraph tree or data.frame with from and to columns.
#' @return Named numeric vector with log_likelihood, log_prior, log_posterior
#'   and posterior_probability. Trees outside hard support receive zero mass.
#' @export
score_tree <- function(fit, tree) {
  if (!inherits(fit, "tconstrain_posterior")) stop("'fit' must come from fit_tree_posterior().")
  if (inherits(tree, "igraph")) {
    if (igraph::is_directed(tree) || !setequal(igraph::V(tree)$name, fit$cluster_ids)) {
      stop("Candidate graph must be undirected with exactly the fitted vertices.")
    }
    tree <- igraph::as_data_frame(tree, what = "edges")
  }
  if (!is.data.frame(tree) || !all(c("from", "to") %in% names(tree))) stop("Tree must have from and to columns.")
  rules <- .tc_rules(data.frame(from = tree$from, to = tree$to, type = rep("must_link", nrow(tree))), fit$cluster_ids)
  if (nrow(tree) != length(fit$cluster_ids) - 1L || nrow(rules) != nrow(tree)) stop("Candidate must be a spanning tree.")
  .tc_support(length(fit$cluster_ids), fit$edges, rules)
  idx <- .tc_rule_indices(rules, fit$edges)
  supported <- all(which(fit$edges$required) %in% idx) && !any(fit$edges$forbidden[idx])
  ll <- sum(fit$edges$log_likelihood[idx])
  lp <- if (supported) sum(fit$edges$log_prior_weight[idx]) - fit$log_prior_normalizer else -Inf
  post <- if (supported) sum(fit$edges$log_posterior_weight[idx]) - fit$log_posterior_normalizer else -Inf
  c(log_likelihood = ll, log_prior = lp, log_posterior = post,
    posterior_probability = exp(post))
}
