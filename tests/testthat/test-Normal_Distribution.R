test_that("Testing normal distribution", {
  library(latentClass)
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
