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
  
  model$add_normal("x1",
                   x,
                   0,
                   1,
                   FALSE)
  testthat::expect_equal(model$get_parameters(),
                         list(x1 = matrix(c(0,1),
                                          nrow = 2,
                                          byrow = TRUE,
                                          dimnames = list(c("mean", "sd"), paste0("class_", 1)))))
  
  testthat::expect_equal(model$log_likelihood(),
                         sum(dnorm(x = x, mean = 0, sd = 1, log = TRUE)))
  
})



test_that("Testing LCA - multi-class normal", {
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
  
  model$add_normal("x1",
                   x,
                   c(0,1,2,3),
                   c(1,1,1,1),
                   TRUE)
  model$get_parameters()
  
  testthat::expect_equal(model$get_parameters(),
                         list(x1 = matrix(c(0,1,2,3,
                                            1,1,1,1),
                                          nrow = 2,
                                          byrow = TRUE,
                                          dimnames = list(c("mean", "sd"), paste0("class_", 1:4)))))
  
  model$expectation_maximization()
  
  testthat::expect_equal(unname(sort(model$get_parameters()$x1["mean",])),
                         c(-3, 1, 3, 5),
                         tolerance = .1)
  testthat::expect_equal(unname(model$get_parameters()$x1["sd", 1]),
                         1,
                         tolerance = .1)
  
  pars_out <- model$get_parameters()$x1
  l <- c()
  for(i in 1:length(x)){
    l_i <- 0
    for(cl in 1:length(model$get_class_probabilities())){
      l_i <- l_i + (model$get_class_probabilities()[cl] * 
                      dnorm(x[i], 
                            pars_out["mean", cl], 
                            pars_out["sd", cl]))
    }
    l <- c(l, l_i)
  }
  testthat::expect_equal(
    model$log_likelihood(),
    sum(log(l)))
  
})
