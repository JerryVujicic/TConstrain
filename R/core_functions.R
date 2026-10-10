#' Construct a hard-constrained MST or a posterior MAP tree
#'
#' Hard mode enforces feasible rules exactly. Probabilistic mode uses the
#' normalized Gaussian tree model in fit_tree_posterior, replacing the historical
#' distance-multiplier heuristic. Time supervision is an independent, symmetric
#' absolute-gap penalty and never orients edges.
#'
#' @param centroids Numeric coordinates with unique cluster row names.
#' @param constraints Data.frame with from, to, type; hard rules when
#'   probabilistic is FALSE and soft rules otherwise.
#' @param time_labels Named numeric times. NA means unknown.
#' @param probabilistic Single nonmissing logical. TRUE returns a posterior MAP tree.
#' @param prior_strength Finite number >= 1; soft tree-prior weight ratio.
#' @param small_eps Deprecated positive finite compatibility argument, ignored.
#' @param large_value Deprecated positive finite compatibility argument, ignored.
#' @param return_format Either igraph or edgelist. The graph stores audit and
#'   posterior attributes; the edge list stores these as R attributes.
#' @param likelihood_scale Fixed positive Gaussian increment standard deviation.
#' @param time_strength Nonnegative absolute time-gap penalty; default 0.
#' @param time_scale Positive divisor for time units.
#' @param hard_constraints Additional exact constraints, also available in soft mode.
#' @return An undirected igraph tree or data.frame with from, to and weight.
#'   The weight is an optimization cost, not pseudotime. In probabilistic mode
#'   it is the negative edge log posterior weight and may be negative.
#'   In hard mode it is Euclidean distance plus the temporal penalty.
#'   Infeasible hard rules produce an error, never an artificial bridge.
#' @export
get_constrained_mst <- function(centroids, constraints = NULL, time_labels = NULL,
                                probabilistic = FALSE, prior_strength = 10,
                                small_eps = 1e-8, large_value = 1e12,
                                return_format = c("igraph", "edgelist"),
                                likelihood_scale = 1, time_strength = 0,
                                time_scale = 1, hard_constraints = NULL) {
  return_format <- match.arg(return_format)
  if (!is.logical(probabilistic) || length(probabilistic) != 1L || is.na(probabilistic)) {
    stop("'probabilistic' must be a single nonmissing logical.")
  }
  .tc_scalar(small_eps, "small_eps")
  .tc_scalar(large_value, "large_value")
  if (!missing(small_eps) || !missing(large_value)) {
    warning("small_eps and large_value are deprecated and ignored; hard rules are exact.", call. = FALSE)
  }
  .tc_scalar(prior_strength, "prior_strength", 1, FALSE)
  .tc_scalar(likelihood_scale, "likelihood_scale")
  .tc_scalar(time_strength, "time_strength", 0, FALSE)
  .tc_scalar(time_scale, "time_scale")
  x <- .tc_coordinates(centroids)
  ids <- rownames(x)
  if (probabilistic) {
    posterior <- fit_tree_posterior(x, constraints, hard_constraints, time_labels,
                                    prior_strength, likelihood_scale, time_strength, time_scale)
    tree <- posterior$map
    edges <- posterior$edges
    cost <- -edges$log_posterior_weight
  } else {
    posterior <- NULL
    # Validate both tables before combining them, including zero-row tables.
    rules <- rbind(.tc_rules(constraints, ids), .tc_rules(hard_constraints, ids))
    rules <- .tc_rules(data.frame(from = ids[rules$i], to = ids[rules$j], type = rules$type), ids)
    edges <- .tc_edges(x, time_labels, time_scale)
    support <- .tc_support(nrow(x), edges, rules)
    cost <- edges$distance + time_strength * edges$time_gap
    selected <- .tc_map(nrow(x), edges, support, cost)
    edges$required <- support$required
    edges$forbidden <- support$forbidden
    edges$eligible <- support$eligible
    edges$in_map <- seq_len(nrow(edges)) %in% selected
    edges$time_penalty <- time_strength * edges$time_gap
    edges$cost <- cost
    tree <- .tc_graph(ids, edges, selected, cost)
  }
  audit <- list(edges = edges, mode = if (probabilistic) "posterior_map" else "hard_mst",
                connected = TRUE, is_tree = TRUE,
                hard_constraints_satisfied = all(edges$in_map[edges$required]) &&
                  !any(edges$in_map[edges$forbidden]))
  # Infinity here is audit metadata for a forbidden edge, never an MST weight.
  cost[edges$forbidden] <- Inf
  igraph::graph_attr(tree, "modified_distance_matrix") <- .tc_matrix(ids, edges, cost)
  igraph::graph_attr(tree, "audit") <- audit
  if (!is.null(posterior)) igraph::graph_attr(tree, "posterior") <- posterior
  if (return_format == "igraph") return(tree)
  out <- igraph::as_data_frame(tree, what = "edges")[, c("from", "to", "weight"), drop = FALSE]
  attr(out, "audit") <- audit
  attr(out, "posterior") <- posterior
  attr(out, "cluster_ids") <- ids
  out
}

#' Construct a constrained tree from a data frame of centroids
#'
#' @param data Data.frame containing an ID column and numeric coordinates only.
#' @param id_col Name of the cluster ID column.
#' @inheritParams get_constrained_mst
#' @return An igraph tree or edge data.frame, as in get_constrained_mst.
#' @export
get_constrained_mst_tidy <- function(data, id_col = "id", constraints = NULL,
                                     time_labels = NULL, probabilistic = FALSE,
                                     prior_strength = 10, small_eps = 1e-8,
                                     large_value = 1e12, return_format = c("igraph", "edgelist"),
                                     likelihood_scale = 1, time_strength = 0,
                                     time_scale = 1, hard_constraints = NULL) {
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!is.character(id_col) || length(id_col) != 1L || is.na(id_col) || !id_col %in% names(data)) stop("Invalid id_col.")
  cent <- as.data.frame(data)
  ids <- as.character(cent[[id_col]])
  .tc_ids(ids, nrow(cent))
  rownames(cent) <- ids
  cent[[id_col]] <- NULL
  args <- list(centroids = cent, constraints = constraints, time_labels = time_labels,
               probabilistic = probabilistic, prior_strength = prior_strength,
               return_format = match.arg(return_format), likelihood_scale = likelihood_scale,
               time_strength = time_strength, time_scale = time_scale, hard_constraints = hard_constraints)
  if (!missing(small_eps)) args$small_eps <- small_eps
  if (!missing(large_value)) args$large_value <- large_value
  do.call(get_constrained_mst, args)
}

#' Run deterministic package smoke checks
#'
#' @param verbose Logical; whether to print a success message.
#' @return TRUE invisibly on success. Does not change the random number state.
#' @export
run_constrained_mst_checks <- function(verbose = TRUE) {
  x <- rbind(C1 = c(0, 0), C2 = c(1, 0), C3 = c(0, 1))
  rule <- data.frame(from = "C1", to = "C2", type = "cannot_link")
  tree <- get_constrained_mst(x, rule)
  stopifnot(!igraph::are_adjacent(tree, "C1", "C2"), igraph::is_tree(tree))
  fit <- fit_tree_posterior(x)
  stopifnot(abs(sum(fit$edges$posterior_probability) - 2) < 1e-10,
            abs(score_tree(fit, fit$map)["posterior_probability"] - fit$map_probability) < 1e-10)
  if (verbose) message("All smoke checks passed. Run the testthat suite for full mathematical checks.")
  invisible(TRUE)
}
