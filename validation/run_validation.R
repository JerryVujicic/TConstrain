# Run from the repository root after installing TConstrain 0.1.0.
library(TConstrain)
source("tests/testthat/helper-enumeration.R")
dir.create("validation/results", showWarnings = FALSE, recursive = TRUE)
set.seed(20261010)

# Exhaustive comparisons: 16 trees for K=4, 125 trees for K=5.
errors <- list()
for (n in 2:5) {
  x <- matrix(rnorm(n * 3), n, 3, dimnames = list(LETTERS[seq_len(n)], NULL))
  fit <- fit_tree_posterior(x, constraints = rule("A", "B"), prior_strength = 7,
                            likelihood_scale = 2,
                            time_labels = setNames(seq_len(n), rownames(x)),
                            time_strength = 0.3, time_scale = 2)
  oracle <- enumerate_trees(fit)
  errors[[n - 1]] <- data.frame(clusters = n, trees = length(oracle$sets),
    prior_z_error = abs(fit$log_prior_normalizer - oracle$log_prior_z),
    posterior_z_error = abs(fit$log_posterior_normalizer - oracle$log_posterior_z),
    marginal_error = max(abs(fit$edges$posterior_probability - oracle$marginal)),
    map_error = abs(fit$map_probability - max(oracle$posterior)))
}
write.csv(do.call(rbind, errors), "validation/results/enumeration.csv", row.names = FALSE)

# Original rebound geometry; reproduce the old heuristic rather than silently
# comparing it to a different objective under the same name.
x <- rbind(Day0 = c(0, 0), Day3 = c(5, 5), Day5 = c(5, 0), Day10 = c(1, 1))
cons <- rbind(rule("Day0", "Day3"), rule("Day0", "Day10", "cannot_link"))
times <- c(Day0 = 0, Day3 = 3, Day5 = 5, Day10 = 10)
base <- fit_tree_posterior(x, likelihood_scale = 5)
sets <- enumerate_trees(base)$sets
old_cost <- base$edges$distance * 10.5
old_cost[base$edges$from == "Day0" & base$edges$to == "Day3"] <- old_cost[base$edges$from == "Day0" & base$edges$to == "Day3"] / 20
old_cost[base$edges$from == "Day0" & base$edges$to == "Day10"] <- old_cost[base$edges$from == "Day0" & base$edges$to == "Day10"] * 20
costs <- vapply(sets, function(k) sum(old_cost[k]), numeric(1))
edge_names <- function(k, edges = base$edges) paste(paste(edges$from[k], edges$to[k], sep = "--"), collapse = "; ")
legacy <- data.frame(edges = vapply(sets, edge_names, character(1)), old_cost = costs,
                     optimal = abs(costs - min(costs)) < 1e-10)
stopifnot(sum(legacy$optimal) == 2L)
write.csv(legacy, "validation/results/rebound_legacy.csv", row.names = FALSE)
models <- list(
  geometry = fit_tree_posterior(x, likelihood_scale = 5),
  soft = fit_tree_posterior(x, cons, likelihood_scale = 5, prior_strength = 20),
  soft_time = fit_tree_posterior(x, cons, time_labels = times, likelihood_scale = 5,
                                prior_strength = 20, time_strength = 1, time_scale = 5))
rebound <- do.call(rbind, lapply(names(models), function(name) {
  fit <- models[[name]]
  data.frame(model = name, fit$edges, map_probability = fit$map_probability,
             log_evidence = fit$log_evidence, entropy = fit$entropy)
}))
write.csv(rebound, "validation/results/rebound_posterior.csv", row.names = FALSE)

grid <- expand.grid(likelihood_scale = c(2, 5, 10), prior_strength = c(1, 5, 20, 100),
                     time_strength = c(0, 1))
for (j in seq_len(nrow(grid))) {
  fit <- fit_tree_posterior(x, cons, time_labels = times, likelihood_scale = grid$likelihood_scale[j],
                            prior_strength = grid$prior_strength[j], time_strength = grid$time_strength[j], time_scale = 5)
  grid$map_probability[j] <- fit$map_probability
  grid$shortcut_probability[j] <- fit$edges$posterior_probability[fit$edges$from == "Day0" & fit$edges$to == "Day10"]
  grid$map_edges[j] <- edge_names(which(fit$edges$in_map), fit$edges)
}
write.csv(grid, "validation/results/sensitivity.csv", row.names = FALSE)

# Simulation-based calibration under the stated generative model, not a
# biological accuracy benchmark: draw a uniform tree, anchor A at zero, draw
# independent N(0,1) increments in two dimensions. No tree edges are supplied.
reference <- fit_tree_posterior(matrix(0, 4, 2, dimnames = list(LETTERS[1:4], NULL)))
trees <- enumerate_trees(reference)$sets
simulation <- vector("list", 400)
coverage <- logical(400)
recovery <- numeric(400)
for (replicate in seq_len(400)) {
  true_index <- sample.int(length(trees), 1)
  true_edges <- trees[[true_index]]
  graph <- igraph::graph_from_data_frame(reference$edges[true_edges, c("from", "to")],
                                        directed = FALSE, vertices = data.frame(name = LETTERS[1:4]))
  x <- matrix(NA_real_, 4, 2, dimnames = list(LETTERS[1:4], NULL))
  x[1, ] <- 0
  visit_order <- as.integer(igraph::bfs(graph, root = "A")$order)
  for (j in visit_order[-1]) {
    parent <- as.integer(igraph::neighbors(graph, j))
    parent <- parent[!is.na(x[parent, 1])]
    stopifnot(length(parent) == 1)
    x[j, ] <- x[parent, ] + rnorm(2)
  }
  fit <- fit_tree_posterior(x, likelihood_scale = 1, prior_strength = 1)
  truth <- as.integer(seq_len(6) %in% true_edges)
  simulation[[replicate]] <- data.frame(replicate = replicate, edge = seq_len(6),
    probability = fit$edges$posterior_probability, truth = truth,
    brier = (fit$edges$posterior_probability - truth)^2)
  post <- vapply(trees, function(k) score_tree(fit, fit$edges[k, , drop = FALSE])["posterior_probability"], numeric(1))
  rank <- order(post, decreasing = TRUE)
  count <- which(cumsum(post[rank]) >= 0.9)[1]
  coverage[replicate] <- true_index %in% rank[seq_len(count)]
  recovery[replicate] <- mean(fit$edges$in_map[true_edges])
}
sim <- do.call(rbind, simulation)
sim$bin <- cut(sim$probability, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
calibration <- do.call(rbind, lapply(split(sim, sim$bin), function(d) {
  data.frame(bin = as.character(d$bin[1]), n_edges = nrow(d),
             mean_probability = mean(d$probability), observed_fraction = mean(d$truth))
}))
write.csv(calibration, "validation/results/calibration.csv", row.names = FALSE)
write.csv(sim, "validation/results/simulation_edges.csv", row.names = FALSE)
summary <- data.frame(replicates = 400, dimensions = 2, clusters = 4,
  mean_edge_brier = mean(sim$brier), uniform_prior_brier = 0.25,
  credible_set_90_coverage = mean(coverage), mean_map_true_edge_recall = mean(recovery))
write.csv(summary, "validation/results/simulation_summary.csv", row.names = FALSE)

# Indicative runtime for this machine; complete dense graphs, not a formal benchmark.
timings <- do.call(rbind, lapply(c(10, 50, 100), function(n) {
  x <- matrix(rnorm(n * 3), n, 3, dimnames = list(paste0("C", seq_len(n)), NULL))
  elapsed <- system.time(fit_tree_posterior(x, likelihood_scale = 3))["elapsed"]
  data.frame(clusters = n, candidate_edges = n * (n - 1) / 2, seconds = unname(elapsed))
}))
write.csv(timings, "validation/results/runtime.csv", row.names = FALSE)
session <- gsub(dirname(dirname(R.home())), "<R environment>", capture.output(sessionInfo()), fixed = TRUE)
writeLines(session, "validation/results/sessionInfo.txt")
print(do.call(rbind, errors))
print(summary)
print(timings)

# Deliberately incorrect prior: a transparency check, not a biological benchmark.
x <- rbind(A = 0, B = 1, C = 2, D = 3)
controls <- list(baseline = fit_tree_posterior(x),
                 correct_rule = fit_tree_posterior(x, rule("A", "B"), prior_strength = 100),
                 incorrect_rule = fit_tree_posterior(x, rule("A", "D"), prior_strength = 10000))
negative <- do.call(rbind, lapply(names(controls), function(name) {
  fit <- controls[[name]]
  data.frame(model = name, long_edge_probability = fit$edges$posterior_probability[fit$edges$from == "A" & fit$edges$to == "D"],
             map_probability = fit$map_probability,
             chain_edges_in_map = sum(fit$edges$in_map & fit$edges$distance == 1))
}))
write.csv(negative, "validation/results/negative_controls.csv", row.names = FALSE)
print(negative)
