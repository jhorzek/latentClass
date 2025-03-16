test_that("testing weighted estimation - normals", {
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

  weights <- runif(n = nrow(normal_data), min = .1, max = 1)

  fit <- latentClass(
    data = normal_data,
    normal = normal(items = c("norm_1", "norm_2", "norm_3")),
    n_classes = 3,
    sample_weights = weights
  )

  # recreate weighted log-likelihood
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
    ind_ll[i] <- weights[i] * log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )

  # mclust also offers a weighted estimation, allowing us to check the results agains mclust
  library(mclust)
  z <- mclust::unmap(sample(1:3, size = nrow(normal_data), replace = TRUE))
  mclust_fit <- me.weighted(
    data = normal_data,
    modelName = "EEI",
    weights = weights,
    z = z
  )

  testthat::expect_equal(
    fit$fit$`log-Likelihood` / mean(weights),
    mclust_fit$loglik,
    tolerance = .01
  )

  testthat::expect_equal(
    fit$fit$BIC,
    -mclust_fit$bic,
    tolerance = .01
  )
})

test_that("testing weighted estimation - categoricals", {
  library(latentClass)
  set.seed(123)
  # only categorical items
  categorical_data <- data.frame(
    cat_1 = factor(c(
      sample(1:3, 200, replace = TRUE, prob = c(.8, .1, .1)),
      sample(1:3, 200, replace = TRUE, prob = c(.3, .3, .4)),
      sample(1:3, 100, replace = TRUE, prob = c(.2, .6, .2))
    )),
    cat_2 = factor(c(
      sample(1:2, 200, replace = TRUE, prob = c(.3, .7)),
      sample(1:2, 200, replace = TRUE, prob = c(.6, .4)),
      sample(1:2, 100, replace = TRUE, prob = c(.1, .9))
    ))
  )

  weights <- runif(n = nrow(categorical_data), min = .1, max = 1)

  fit <- latentClass(
    data = categorical_data,
    categorical = categorical(items = c("cat_1", "cat_2")),
    n_classes = 3,
    sample_weights = weights
  )

  # recreate weighted log-likelihood
  ind_ll <- rep(0, nrow(categorical_data))
  for (i in 1:nrow(categorical_data)) {
    ind_L <- fit$class_probabilities
    for (cl in 1:length(ind_L)) {
      for (j in 1:ncol(categorical_data)) {
        ind_L[cl] <- ind_L[cl] *
          fit$estimates[[colnames(categorical_data)[j]]][
            categorical_data[i, j],
            cl
          ]
      }
    }
    ind_ll[i] <- weights[i] * log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )
})

test_that("testing weighted estimation - normals and categoricals", {
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
  # only categorical items
  categorical_data <- data.frame(
    cat_1 = factor(c(
      sample(1:3, 200, replace = TRUE, prob = c(.8, .1, .1)),
      sample(1:3, 200, replace = TRUE, prob = c(.3, .3, .4)),
      sample(1:3, 100, replace = TRUE, prob = c(.2, .6, .2))
    )),
    cat_2 = factor(c(
      sample(1:2, 200, replace = TRUE, prob = c(.3, .7)),
      sample(1:2, 200, replace = TRUE, prob = c(.6, .4)),
      sample(1:2, 100, replace = TRUE, prob = c(.1, .9))
    ))
  )

  data <- cbind(normal_data, categorical_data)

  weights <- runif(n = nrow(normal_data), min = .1, max = 1)

  fit <- latentClass(
    data = data,
    categorical = categorical(items = c("cat_1", "cat_2")),
    normal = normal(items = c("norm_1", "norm_2", "norm_3")),
    n_classes = 3,
    sample_weights = weights
  )

  # recreate weighted log-likelihood
  ind_ll <- rep(0, nrow(data))
  for (i in 1:nrow(data)) {
    ind_L <- fit$class_probabilities
    for (cl in 1:length(ind_L)) {
      for (j in 1:ncol(data)) {
        if (is.factor(data[[j]])) {
          ind_L[cl] <- ind_L[cl] *
            fit$estimates[[colnames(data)[j]]][data[i, j], cl]
        } else {
          ind_L[cl] <- ind_L[cl] *
            dnorm(
              data[i, j],
              mean = fit$estimates[[colnames(data)[j]]]["mean", cl],
              sd = fit$estimates[[colnames(data)[j]]]["sd", cl]
            )
        }
      }
    }
    ind_ll[i] <- weights[i] * log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )
})
