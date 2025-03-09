test_that("Testing LCA - single class normal", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- rnorm(100)
  
  model <- latentClass:::LCMR$new(1,
                                  length(x))
  testthat::expect_equal(model$get_n_classes(), 1)
  testthat::expect_equal(model$get_class_probabilities(), 1)
  testthat::expect_equal(model$get_n_persons(), 100)
  
  model$add_normal(1,
                   x,
                   c("mu", "log_sd"),
                   c(0, 0))
  testthat::expect_equal(model$get_parameters(),
                         c("mu" = 0, "log_sd" = 0))
  
  model$update_responsibilities()
  
  testthat::expect_equal(model$get_responsibilities(),
                         matrix(1,
                                nrow = 100,
                                ncol = 1))
  
  pars <- model$get_parameters()
  
  testthat::expect_equal(model$log_likelihood(c("mu", "log_sd"),
                                              c(0, 0)),
                         sum(dnorm(x = x, mean = 0, sd = 1, log = TRUE)))
  
  testthat::expect_equal(
    unname(model$gradients(c("mu", "log_sd"),
                           c(0, 0))),
    numDeriv::grad(func = function(par){
      sum(dnorm(x = x, mean = par[1], sd = exp(par[2]), log = TRUE))
    }, x = c(0, 0)))
  
})

