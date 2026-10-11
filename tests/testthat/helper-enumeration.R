# Independent oracle: enumerate every edge subset, then test graph connectivity.
# This deliberately does not use the package's contraction/MAP algorithms.
enumerate_trees <- function(fit) {
  e <- fit$edges
  n <- length(fit$cluster_ids)
  sets <- if (n == 1L) list(integer()) else combn(seq_len(nrow(e)), n - 1L, simplify = FALSE)
  sets <- Filter(function(k) {
    g <- igraph::graph_from_data_frame(e[k, c("from", "to"), drop = FALSE],
                                      directed = FALSE, vertices = data.frame(name = fit$cluster_ids))
    igraph::is_connected(g) && all(which(e$required) %in% k) && !any(e$forbidden[k])
  }, sets)
  lp <- vapply(sets, function(k) sum(e$log_prior_weight[k]), numeric(1))
  ll <- vapply(sets, function(k) sum(e$log_likelihood[k]), numeric(1))
  lse <- function(v) max(v) + log(sum(exp(v - max(v))))
  z0 <- lse(lp)
  z1 <- lse(lp + ll)
  pp <- exp(lp + ll - z1)
  p0 <- exp(lp - z0)
  marginal <- prior <- numeric(nrow(e))
  for (k in seq_along(sets)) {
    marginal[sets[[k]]] <- marginal[sets[[k]]] + pp[k]
    prior[sets[[k]]] <- prior[sets[[k]]] + p0[k]
  }
  list(sets = sets, log_prior_z = z0, log_posterior_z = z1, posterior = pp,
       prior = p0, marginal = marginal, prior_marginal = prior,
       entropy = -sum(pp * log(pp)))
}

rule <- function(a, b, type = "must_link") data.frame(from = a, to = b, type = type)
