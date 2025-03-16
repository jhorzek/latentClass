test_that("Testing LCA - single class normal", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- rnorm(100)

  model <- latentClass:::LCMR$new(1, length(x))
  testthat::expect_equal(model$get_n_classes(), 1)
  testthat::expect_equal(model$get_class_probabilities(), 1)
  testthat::expect_equal(model$get_n_persons(), 100)

  model$add_normal("x1", x, 0, 1, FALSE)
  testthat::expect_equal(
    model$get_parameters(),
    list(
      x1 = matrix(
        c(0, 1),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("mean", "sd"), paste0("class_", 1))
      )
    )
  )

  testthat::expect_equal(
    model$log_likelihood(),
    sum(dnorm(x = x, mean = 0, sd = 1, log = TRUE))
  )
})


test_that("Testing LCA - multi-class normal", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- c(
    rnorm(500, -3, exp(0)),
    rnorm(500, 3, exp(0)),
    rnorm(500, 5, exp(0)),
    rnorm(500, 1, exp(0))
  )
  model <- latentClass:::LCMR$new(rep(1, 4) / 4, length(x))
  testthat::expect_equal(model$get_n_classes(), 4)
  testthat::expect_equal(model$get_class_probabilities(), rep(1, 4) / 4)
  testthat::expect_equal(model$get_n_persons(), 2000)

  model$add_normal("x1", x, c(0, 1, 2, 3), c(1, 1, 1, 1), TRUE)
  model$get_parameters()

  testthat::expect_equal(
    model$get_parameters(),
    list(
      x1 = matrix(
        c(0, 1, 2, 3, 1, 1, 1, 1),
        nrow = 2,
        byrow = TRUE,
        dimnames = list(c("mean", "sd"), paste0("class_", 1:4))
      )
    )
  )

  model$expectation_maximization(1000, 1e-10)

  testthat::expect_equal(
    unname(sort(model$get_parameters()$x1["mean", ])),
    c(-3, 1, 3, 5),
    tolerance = .1
  )
  testthat::expect_equal(
    unname(model$get_parameters()$x1["sd", 1]),
    1,
    tolerance = .1
  )

  pars_out <- model$get_parameters()$x1
  l <- c()
  for (i in 1:length(x)) {
    l_i <- 0
    for (cl in 1:length(model$get_class_probabilities())) {
      l_i <- l_i +
        (model$get_class_probabilities()[cl] *
          dnorm(x[i], pars_out["mean", cl], pars_out["sd", cl]))
    }
    l <- c(l, l_i)
  }
  testthat::expect_equal(
    model$log_likelihood(),
    sum(log(l))
  )
})


test_that("testing unweighted estimation - normals", {
  library(latentClass)
  set.seed(123)
  # only gaussian items
  normal_data <- data.frame(
    # no difference in first item
    norm_1 = c(
      rnorm(n = 200, mean = 0, sd = 1),
      rnorm(n = 200, mean = 0, sd = 1),
      rnorm(n = 100, mean = 0, sd = 1)
    ),
    # difference in second item
    norm_2 = c(
      rnorm(n = 200, mean = 3, sd = 1),
      rnorm(n = 200, mean = -3, sd = 1),
      rnorm(n = 100, mean = 0, sd = 1)
    ),
    # difference in third item
    norm_3 = c(
      rnorm(n = 200, mean = .5, sd = .3),
      rnorm(n = 200, mean = .1, sd = .3),
      rnorm(n = 100, mean = .2, sd = .3)
    )
  )

  fit <- latentClass(
    data = normal_data,
    normal = normal(items = c("norm_1", "norm_2", "norm_3")),
    n_classes = 3
  )

  # recreate log-likelihood
  ind_ll <- rep(0, nrow(normal_data))
  for (i in 1:nrow(normal_data)) {
    ind_L <- fit$class_probabilities
    for (cl in 1:length(ind_L)) {
      for (j in 1:ncol(normal_data)) {
        ind_L[cl] <- ind_L[cl] *
          dnorm(
            normal_data[i, j],
            mean = fit$estimates[[colnames(normal_data)[j]]]["mean", cl],
            sd = fit$estimates[[colnames(normal_data)[j]]]["sd", cl]
          )
      }
    }
    ind_ll[i] <- log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )

  # Check the results agains mclust
  library(mclust)
  z <- mclust::unmap(sample(1:3, size = nrow(normal_data), replace = TRUE))
  mclust_fit <- me.weighted(
    data = normal_data,
    modelName = "EEI",
    weights = rep(1, nrow(normal_data)),
    z = z
  )

  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    mclust_fit$loglik,
    tolerance = .01
  )

  testthat::expect_equal(
    fit$fit$BIC,
    -mclust_fit$bic,
    tolerance = .01
  )
})


test_that("testing unweighted estimation - normals with free sd", {
  library(latentClass)
  set.seed(123)
  # only gaussian items
  normal_data <- data.frame(
    # no difference in first item
    norm_1 = c(
      rnorm(n = 200, mean = 0, sd = 1),
      rnorm(n = 200, mean = 0, sd = 1),
      rnorm(n = 100, mean = 0, sd = 1)
    ),
    # difference in second item
    norm_2 = c(
      rnorm(n = 200, mean = 3, sd = 1),
      rnorm(n = 200, mean = -3, sd = 1),
      rnorm(n = 100, mean = 0, sd = 1)
    ),
    # difference in third item
    norm_3 = c(
      rnorm(n = 200, mean = .5, sd = 1.3),
      rnorm(n = 200, mean = .1, sd = .5),
      rnorm(n = 100, mean = .2, sd = .2)
    )
  )

  fit <- latentClass(
    data = normal_data,
    normal = normal(items = c("norm_1", "norm_2", "norm_3"), sd_equal = FALSE),
    n_classes = 3
  )

  # recreate log-likelihood
  ind_ll <- rep(0, nrow(normal_data))
  for (i in 1:nrow(normal_data)) {
    ind_L <- fit$class_probabilities
    for (cl in 1:length(ind_L)) {
      for (j in 1:ncol(normal_data)) {
        ind_L[cl] <- ind_L[cl] *
          dnorm(
            normal_data[i, j],
            mean = fit$estimates[[colnames(normal_data)[j]]]["mean", cl],
            sd = fit$estimates[[colnames(normal_data)[j]]]["sd", cl]
          )
      }
    }
    ind_ll[i] <- log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )

  # Check the results agains mclust
  library(mclust)
  z <- mclust::unmap(sample(1:3, size = nrow(normal_data), replace = TRUE))
  mclust_fit <- me.weighted(
    data = normal_data,
    modelName = "VVI",
    weights = rep(1, nrow(normal_data)),
    z = z
  )

  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    mclust_fit$loglik,
    tolerance = .01
  )

  testthat::expect_equal(
    fit$fit$BIC,
    -mclust_fit$bic,
    tolerance = .01
  )

  stop("Parameter estimates differ between mclust and latentClass.")
})
