test_that("Testing normal distribution", {
  library(latentClass)
  set.seed(234)
  x <- rnorm(100)
  n <- latentClass:::NormalDistribution$new(x,
                                            c("mu", "l_sigma"),
                                            c(1,0),
                                            c(TRUE, TRUE))
  
  testthat::expect_equal(n$get_log_likelihood(c("mu", "l_sigma"),
                                              c(1,0)),
                         sum(dnorm(x = x, mean = 1, sd = 1, log = TRUE)))
  
  
  testthat::expect_equal(n$get_individual_log_likelihood(c("mu", "l_sigma"),
                                                         c(1,0)),
                         dnorm(x = x, mean = 1, sd = 1, log = TRUE))
  
  testthat::expect_equal(
    n$get_gradients(c("mu", "l_sigma"),
                    c(1,0)),
    numDeriv::grad(func = function(par){
      sum(dnorm(x = x, mean = par[1], sd = exp(par[2]), log = TRUE))
    }, x = c(1,0))
  )
  
})


test_that("Testing normal distribution with weights", {
  library(latentClass)
  set.seed(234)
  x <- rnorm(100)
  w <- runif(n = length(x), min = .1, max = 3)
  n <- latentClass:::NormalDistribution$new(x,
                                            c("mu", "l_sigma"),
                                            c(1,0),
                                            c(TRUE, TRUE))
  
  testthat::expect_equal(n$get_log_likelihood_w(c("mu", "l_sigma"),
                                                c(1,0),
                                                w),
                         sum(w*dnorm(x = x, mean = 1, sd = 1, log = TRUE)))
  
  
  testthat::expect_equal(n$get_individual_log_likelihood_w(c("mu", "l_sigma"),
                                                           c(1,0),
                                                           w),
                         w*dnorm(x = x, mean = 1, sd = 1, log = TRUE))
  
  testthat::expect_equal(
    n$get_gradients_w(c("mu", "l_sigma"),
                      c(1,0),
                      w),
    numDeriv::grad(func = function(par){
      sum(w*dnorm(x = x, mean = par[1], sd = exp(par[2]), log = TRUE))
    }, x = c(1,0))
  )
  
})
