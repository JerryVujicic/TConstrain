test_that("normalizers, marginals, MAP, evidence and tree scoring match enumeration", {
  for (n in 2:5) {
    x <- cbind(seq_len(n) / 3, sin(seq_len(n)))
    rownames(x) <- LETTERS[seq_len(n)]
    for (hard in list(NULL, rule("A", "B"), rule("A", "B", "cannot_link"))) {
      if (n == 2L && !is.null(hard) && hard$type == "cannot_link") next
      fit <- fit_tree_posterior(x, rule("A", "B"), hard_constraints = hard,
                                time_labels = setNames(seq_len(n), rownames(x)),
                                prior_strength = 3, likelihood_scale = 2,
                                time_strength = 0.2, time_scale = 2)
      oracle <- enumerate_trees(fit)
      expect_equal(fit$log_prior_normalizer, oracle$log_prior_z, tolerance = 1e-10)
      expect_equal(fit$log_posterior_normalizer, oracle$log_posterior_z, tolerance = 1e-10)
      expect_equal(fit$log_evidence, oracle$log_posterior_z - oracle$log_prior_z, tolerance = 1e-10)
      expect_equal(fit$edges$posterior_probability, oracle$marginal, tolerance = 1e-10)
      expect_equal(fit$edges$prior_probability, oracle$prior_marginal, tolerance = 1e-10)
      expect_equal(sum(fit$edges$posterior_probability), n - 1, tolerance = 1e-10)
      expect_equal(fit$map_probability, max(oracle$posterior), tolerance = 1e-10)
      expect_equal(fit$entropy, oracle$entropy, tolerance = 1e-10)
      expect_equal(sum(oracle$posterior), 1)
      expect_equal(sum(oracle$prior), 1)
      scored <- vapply(oracle$sets, function(k) score_tree(fit, fit$edges[k, , drop = FALSE])["posterior_probability"], numeric(1))
      expect_equal(scored, oracle$posterior, ignore_attr = TRUE)
    }
  }
})

test_that("uniform coincident geometry follows Cayley's formula", {
  x <- matrix(0, 5, 2, dimnames = list(LETTERS[1:5], NULL))
  fit <- fit_tree_posterior(x, prior_strength = 1)
  expect_equal(fit$log_prior_normalizer, log(5^3))
  expect_equal(fit$map_probability, 1 / 125)
  expect_equal(fit$edges$posterior_probability, rep(2 / 5, 10))
  expect_equal(fit$entropy, log(125))
})

test_that("Gaussian constants and tree-prior exchange ratios are correct", {
  x <- rbind(A = c(0, 0), B = c(1, 0), C = c(0, 1))
  fit <- fit_tree_posterior(x, rule("A", "B"), prior_strength = 4, likelihood_scale = 2)
  expect_equal(fit$edges$log_likelihood, -log(2 * pi * 4) - fit$edges$distance^2 / 8)
  a <- score_tree(fit, rule(c("A", "B"), c("B", "C")))
  b <- score_tree(fit, rule(c("A", "B"), c("C", "C")))
  expect_equal(unname(a["log_prior"] - b["log_prior"]), log(4))
  expect_equal(unname(a["log_posterior"] - b["log_posterior"]),
               unname(a["log_prior"] + a["log_likelihood"] - b["log_prior"] - b["log_likelihood"]))
})

test_that("required contraction preserves parallel-edge probabilities", {
  x <- rbind(A = 0, B = 0, C = 0, D = 0)
  fit <- fit_tree_posterior(x, hard_constraints = rule(c("A", "C"), c("B", "D")))
  expect_equal(fit$edges$posterior_probability[fit$edges$required], c(1, 1))
  expect_equal(fit$edges$posterior_probability[!fit$edges$required], rep(0.25, 4))
  expect_equal(fit$map_probability, 0.25)
  forced <- fit_tree_posterior(x, hard_constraints = rule(c("A", "B", "C"), c("B", "C", "D")))
  expect_equal(forced$map_probability, 1)
  expect_equal(forced$log_prior_normalizer, 0)
  expect_equal(forced$entropy, 0)
})

test_that("scale and translation transformations preserve posterior when specified consistently", {
  x <- rbind(A = c(0, 1), B = c(1, 2), C = c(3, 1))
  a <- fit_tree_posterior(x, likelihood_scale = 2)
  b <- fit_tree_posterior(10 * x + 17, likelihood_scale = 20)
  expect_equal(a$edges$posterior_probability, b$edges$posterior_probability)
  expect_equal(b$log_evidence - a$log_evidence, -4 * log(10))
  c <- fit_tree_posterior(x[c(3, 1, 2), ], likelihood_scale = 2)
  expect_equal(c$log_evidence, a$log_evidence)
})

test_that("posterior includes all tied optimal trees and responds to soft evidence", {
  x <- rbind(A = c(0, 0), B = c(1, 0), C = c(0, 1))
  a <- fit_tree_posterior(x)
  expect_lt(a$map_probability, 0.5)
  b <- fit_tree_posterior(x, rule("A", "B"), prior_strength = 10)
  k <- which(a$edges$from == "A" & a$edges$to == "B")
  expect_gt(b$edges$posterior_probability[k], a$edges$posterior_probability[k])
  expect_lt(b$edges$posterior_probability[k], 1)
  c <- fit_tree_posterior(x, rule("A", "B", "cannot_link"), prior_strength = 10)
  expect_lt(c$edges$posterior_probability[k], a$edges$posterior_probability[k])
})

test_that("unsupported candidate trees have zero posterior mass", {
  x <- rbind(A = 0, B = 1, C = 2)
  fit <- fit_tree_posterior(x, hard_constraints = rule("A", "B", "cannot_link"))
  s <- score_tree(fit, rule(c("A", "B"), c("B", "C")))
  expect_equal(unname(s["posterior_probability"]), 0)
  expect_true(is.finite(s["log_likelihood"]))
  expect_error(score_tree(fit, rule("A", "B")), "spanning tree")
})

test_that("numerically unsupported distributions fail instead of fabricating certainty", {
  x <- rbind(A = 0, B = 1, C = 1000)
  expect_error(fit_tree_posterior(x, likelihood_scale = 0.001), "underflowed|ill-conditioned")
  expect_error(fit_tree_posterior(x, likelihood_scale = 1e-300), "Non-finite")
})
