test_that("testing mixed data", {
  library(latentClass)
  set.seed(123)
  # only gaussian items
  normal_data <- data.frame(
    # no difference in first item
    norm_1 = c(
      rnorm(n = 2000, mean = 0, sd = 1),
      rnorm(n = 2000, mean = 0, sd = 1),
      rnorm(n = 1000, mean = 0, sd = 1)
    ),
    # difference in second item
    norm_2 = c(
      rnorm(n = 2000, mean = 3, sd = 1),
      rnorm(n = 2000, mean = -3, sd = 1),
      rnorm(n = 1000, mean = 0, sd = 1)
    ),
    # difference in third item
    norm_3 = c(
      rnorm(n = 2000, mean = .5, sd = .3),
      rnorm(n = 2000, mean = .1, sd = .3),
      rnorm(n = 1000, mean = .2, sd = .3)
    )
  )
  # only categorical items
  categorical_data <- data.frame(
    cat_1 = factor(c(
      sample(1:3, 2000, replace = TRUE, prob = c(.8, .1, .1)),
      sample(1:3, 2000, replace = TRUE, prob = c(.3, .3, .4)),
      sample(1:3, 1000, replace = TRUE, prob = c(.2, .6, .2))
    )),
    cat_2 = factor(c(
      sample(1:2, 2000, replace = TRUE, prob = c(.3, .7)),
      sample(1:2, 2000, replace = TRUE, prob = c(.6, .4)),
      sample(1:2, 1000, replace = TRUE, prob = c(.1, .9))
    ))
  )

  data <- cbind(normal_data, categorical_data)

  fit <- latentClass(
    data = data,
    categorical = categorical(items = c("cat_1", "cat_2")),
    normal = normal(items = c("norm_1", "norm_2", "norm_3")),
    n_classes = 3
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
    ind_ll[i] <- log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )

  testthat::expect_equal(
    unname(sort(fit$estimates$cat_1[1, ])),
    c(.2, .3, .8),
    tolerance = .1
  )

  testthat::expect_equal(
    unname(sort(fit$estimates$cat_1[2, ])),
    c(.1, .3, .6),
    tolerance = .1
  )

  testthat::expect_equal(
    unname(sort(fit$estimates$cat_1[3, ])),
    c(.1, .2, .4),
    tolerance = .1
  )

  testthat::expect_equal(
    unname(sort(fit$estimates$norm_2["mean", ])),
    c(-3, 0, 3),
    tolerance = .1
  )
})
