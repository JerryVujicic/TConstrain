test_that("hard constraints are literal and infeasibility is explicit", {
  x <- rbind(A = c(0, 0), B = c(1, 0), C = c(0, 2), D = c(1, 2))
  cons <- rbind(rule("A", "D"), rule("A", "B", "cannot_link"))
  g <- get_constrained_mst(x, cons)
  expect_true(igraph::are_adjacent(g, "A", "D"))
  expect_false(igraph::are_adjacent(g, "A", "B"))
  expect_true(igraph::graph_attr(g, "audit")$hard_constraints_satisfied)
  expect_error(get_constrained_mst(x[1:2, ], rule("A", "B", "cannot_link")), "disconnect")
  expect_error(fit_tree_posterior(x, hard_constraints = rule(c("A", "B", "A"), c("B", "C", "C"))), "cycle")
  expect_error(get_constrained_mst(x, rbind(rule("A", "B"), rule("B", "A", "cannot_link"))), "Conflicting")
})

test_that("hard MST agrees with independent enumeration of feasible trees", {
  x <- rbind(A = c(0, 0), B = c(2, 0), C = c(0, 1), D = c(2, 2))
  hard <- rbind(rule("A", "D"), rule("A", "B", "cannot_link"))
  fit <- fit_tree_posterior(x, hard_constraints = hard)
  trees <- enumerate_trees(fit)$sets
  expected <- min(vapply(trees, function(k) sum(fit$edges$distance[k]), numeric(1)))
  g <- get_constrained_mst(x, hard)
  expect_equal(sum(igraph::E(g)$weight), expected)
})

test_that("reverse duplicates do not multiply soft evidence", {
  x <- rbind(A = 0, B = 1, C = 2)
  one <- fit_tree_posterior(x, rule("A", "B"))
  expect_warning(two <- fit_tree_posterior(x, rule(c("A", "B"), c("B", "A"))), "collapsed")
  expect_equal(one$edges, two$edges)
  expect_error(fit_tree_posterior(x, rbind(rule("A", "B"), rule("B", "A", "cannot_link"))), "Conflicting")
})

test_that("zero distance edges, singleton clusters and one dimension are retained", {
  x <- rbind(A = 0, B = 0, C = 1)
  g <- get_constrained_mst(x)
  expect_true(igraph::are_adjacent(g, "A", "B"))
  expect_equal(min(igraph::E(g)$weight), 0)
  expect_equal(dim(get_centroids(matrix(1:4, ncol = 1), c("C2", "C2", "C1", "C1"))), c(2L, 1L))
  one <- infer_trajectory(matrix(c(1, 2), ncol = 1), c("A", "A"), probabilistic = TRUE)
  expect_equal(one$root, "A")
  expect_equal(one$lineages[[1]], "A")
  expect_equal(one$posterior$map_probability, 1)
  expect_equal(one$posterior$log_evidence, 0)
  expect_equal(igraph::ecount(one$mst), 0)
})

test_that("time gaps and pairwise strength are separate and missing times are neutral", {
  x <- rbind(A = 0, B = 1, C = 2)
  times <- c(A = 0, B = 1, C = 3)
  a <- fit_tree_posterior(x, time_labels = times, time_strength = 2, time_scale = 2)
  b <- fit_tree_posterior(x, time_labels = times, time_strength = 2, time_scale = 2, prior_strength = 100)
  expect_equal(a$edges$time_log_prior, c(-1, -3, -2))
  expect_equal(a$edges$posterior_probability, b$edges$posterior_probability)
  c <- fit_tree_posterior(x, time_labels = c(A = 0, B = NA_real_), time_strength = 2)
  expect_equal(c$edges$time_log_prior, c(0, 0, 0))
  expect_equal(fit_tree_posterior(x, time_labels = times)$edges$posterior_probability,
               fit_tree_posterior(x)$edges$posterior_probability)
})

test_that("input errors fail clearly", {
  x <- rbind(A = 0, B = 1, C = 2)
  for (bad in list(NA_real_, Inf, NaN, c(1, 2), -1)) {
    expect_error(fit_tree_posterior(x, likelihood_scale = bad), "likelihood_scale")
  }
  expect_error(get_constrained_mst(x, probabilistic = NA), "nonmissing logical")
  expect_error(fit_tree_posterior(x, time_strength = Inf), "time_strength")
  expect_error(fit_tree_posterior(x, prior_strength = 0.5), "prior_strength")
  expect_error(fit_tree_posterior(x, rule("A", "Z")), "unknown")
  expect_error(fit_tree_posterior(x, rule("A", "A")), "Self")
  expect_error(fit_tree_posterior(x, time_labels = c(A = Inf)), "finite")
  expect_error(fit_tree_posterior(x, time_labels = c(A = 1, A = 2)), "unique")
  expect_error(get_centroids(matrix(c(1, NA), 2), c("A", "B")), "finite")
  expect_error(get_centroids(x, c("A", NA, "C")), "nonmissing")
  expect_error(get_centroids(x, c("A", "", "C")), "nonempty")
  rownames(x) <- c("A", "A", "C")
  expect_error(fit_tree_posterior(x), "unique")
})

test_that("wrappers carry the posterior, support negative costs, and preserve RNG state", {
  x <- rbind(C1 = c(0, 0), C2 = c(0.1, 0), C10 = c(0.2, 0))
  g <- infer_trajectory(x, rownames(x), start_cluster = "C1", probabilistic = TRUE,
                        likelihood_scale = 0.1)
  expect_true(all(igraph::E(g$mst)$weight < 0))
  expect_equal(g$lineages[[1]], c("C1", "C2", "C10"))
  expect_s3_class(g$posterior, "tconstrain_posterior")
  tidy <- get_constrained_mst_tidy(data.frame(id = rownames(x), x), probabilistic = TRUE,
                                   likelihood_scale = 0.1, return_format = "edgelist")
  expect_equal(attr(tidy, "posterior")$map_probability, g$posterior$map_probability)
  expect_error(infer_trajectory(x, rownames(x)), "start_cluster")
  expect_equal(infer_trajectory(x, rownames(x), time_labels = c(C1 = 0, C2 = 1, C10 = 2))$root, "C1")
  set.seed(13)
  before <- .Random.seed
  expect_true(run_constrained_mst_checks(FALSE))
  expect_identical(.Random.seed, before)
})

test_that("hard rules survive arbitrary geometric scales and cannot links cut support", {
  x <- rbind(A = 0, B = 1e20, C = 2e20)
  g <- get_constrained_mst(x, rule("A", "C"))
  expect_true(igraph::are_adjacent(g, "A", "C"))
  cut <- rule(c("A", "A"), c("B", "C"), "cannot_link")
  expect_error(get_constrained_mst(x, cut), "disconnect")
  expect_warning(get_constrained_mst(x, small_eps = 1e-8), "deprecated")
})

test_that("hard restrictions dominate soft preferences in a mixed model", {
  x <- rbind(A = 0, B = 1, C = 2)
  fit <- fit_tree_posterior(x, rule("A", "B"),
                            hard_constraints = rule("A", "B", "cannot_link"), prior_strength = 100)
  expect_equal(fit$edges$posterior_probability[fit$edges$forbidden], 0)
  expect_equal(fit$map_probability, 1)
})
