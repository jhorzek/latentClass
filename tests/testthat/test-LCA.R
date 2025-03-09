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



test_that("Testing LCA - single class normal", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- c(rnorm(500, -3, exp(0)),
         rnorm(500, 3, exp(0)),
         rnorm(500, 5, exp(0)),
         rnorm(500, 1, exp(0)))
  model <- latentClass:::LCMR$new(rep(1, 4)/4,
                                  length(x))
  testthat::expect_equal(model$get_n_classes(),
                         4)
  testthat::expect_equal(model$get_class_probabilities(),
                         rep(1, 4)/4)
  testthat::expect_equal(model$get_n_persons(),
                         2000)
  for(i in 1:4){
    model$add_normal(i,
                     x,
                     c(paste0("mu_", i), paste0("log_sd")),
                     c(i, 0))
  }
  
  testthat::expect_equal(model$get_parameters(),
                         c("mu_1" = 1,
                           "log_sd" = 0,
                           "mu_2" = 2,
                           "mu_3" = 3,
                           "mu_4" = 4))
  latentClass:::expectation_maximization(model = model, 
                                         conv_crit = 1e-10, 
                                         max_iter = 1000, 
                                         use_ad_gradients = TRUE)
  
  testthat::expect_equal(unname(sort(model$get_parameters()[grepl(pattern = "mu_",
                                                                  x = names(model$get_parameters()))])),
                         c(-3, 1, 3, 5),
                         tolerance = .1)
  testthat::expect_equal(unname(model$get_parameters()["log_sd"]),
                         0,
                         tolerance = .1)
  
  pars_out <- model$get_parameters()
  l <- c()
  for(i in 1:length(x)){
    l_i <- 0
    for(cl in 1:length(model$get_class_probabilities())){
      l_i <- l_i + (model$get_class_probabilities()[cl] * 
                      dnorm(x[i], 
                            pars_out[paste0("mu_", cl)], 
                            exp(pars_out["log_sd"])))
    }
    l <- c(l, l_i)
  }
  testthat::expect_equal(
    model$log_likelihood(names(pars_out),
                         pars_out),
    sum(log(l)))
  
})
