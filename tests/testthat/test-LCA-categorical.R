test_that("Testing LCA - single class categorical", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- sample(0:3, size = 1000, replace = TRUE, prob = c(.2, .2, .4, .2))

  model <- latentClass:::LCMR$new(1, length(x))
  testthat::expect_equal(model$get_n_classes(), 1)
  testthat::expect_equal(model$get_class_probabilities(), 1)
  testthat::expect_equal(model$get_n_persons(), 1000)

  model$add_categorical(
    "x1",
    x,
    # starting values for every class
    matrix(c(.1, .2, .3, .4), ncol = 1)
  )
  testthat::expect_equal(
    model$get_parameters(),
    list(
      x1 = matrix(
        c(0.1, 0.2, 0.3, 0.4),
        nrow = 4,
        byrow = TRUE,
        dimnames = list(as.character(0:3), paste0("class_", 1))
      )
    )
  )

  model$expectation_maximization(1000, 1e-10)

  testthat::expect_true(all(
    abs(
      unname(sort(model$get_parameters()$x1[, 1])) -
        unname(sort((table(x) / length(x))))
    ) <
      1e-7
  ))
})


test_that("Testing LCA - multi-class categorical", {
  # For a single class, the expected log-likelihood and the log-likelihood
  # are identical
  library(latentClass)
  set.seed(234)
  x <- c(
    sample(0:3, size = 500, replace = TRUE, prob = c(.2, .2, .4, .2)),
    sample(0:3, size = 500, replace = TRUE, prob = c(.6, .1, .1, .2))
  )
  y <- c(
    sample(0:3, size = 500, replace = TRUE, prob = c(.6, .1, .1, .2)),
    sample(0:3, size = 500, replace = TRUE, prob = c(.2, .2, .4, .2))
  )
  model <- latentClass:::LCMR$new(c(.5, .5), length(x))
  testthat::expect_equal(model$get_n_classes(), 2)
  testthat::expect_equal(model$get_class_probabilities(), c(.5, .5))
  testthat::expect_equal(model$get_n_persons(), 1000)

  model$add_categorical(
    "x1",
    x,
    # starting values for every class
    matrix(c(.1, .2, .3, .4, .2, .2, .3, .3), ncol = 2, byrow = TRUE)
  )
  model$get_parameters()

  testthat::expect_equal(
    model$get_parameters(),
    list(
      x1 = matrix(
        c(.1, .2, .3, .4, .2, .2, .3, .3),
        ncol = 2,
        byrow = TRUE,
        dimnames = list(as.character(0:3), paste0("class_", 1:2))
      )
    )
  )
  model$add_categorical(
    "y1",
    y,
    # starting values for every class
    matrix(c(.2, .2, .3, .3, .1, .2, .3, .4), ncol = 2, byrow = TRUE)
  )
  model$get_parameters()
  model$expectation_maximization(1000, 1e-7)

  # testthat::expect_equal(unname(sort(model$get_parameters()$x1["mean",])),
  #                        c(-3, 1, 3, 5),
  #                        tolerance = .1)
  # testthat::expect_equal(unname(model$get_parameters()$x1["sd", 1]),
  #                        1,
  #                        tolerance = .1)
  #
  # pars_out <- model$get_parameters()$x1
  # l <- c()
  # for(i in 1:length(x)){
  #   l_i <- 0
  #   for(cl in 1:length(model$get_class_probabilities())){
  #     l_i <- l_i + (model$get_class_probabilities()[cl] *
  #                     dnorm(x[i],
  #                           pars_out["mean", cl],
  #                           pars_out["sd", cl]))
  #   }
  #   l <- c(l, l_i)
  # }
  # testthat::expect_equal(
  #   model$log_likelihood(),
  #   sum(log(l)))
})

test_that("testing unweighted estimation - categoricals", {
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

  fit <- latentClass(
    data = categorical_data,
    categorical = categorical(items = c("cat_1", "cat_2")),
    n_classes = 3
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
    ind_ll[i] <- log(sum(ind_L))
  }
  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    sum(ind_ll)
  )

  # Categorical models can also be estimated with poLCA
  library(poLCA)
  fit_polca <- poLCA::poLCA(
    formula = cbind(cat_1, cat_2) ~ 1,
    data = categorical_data,
    nclass = 3
  )

  testthat::expect_equal(
    fit$fit$`log-Likelihood`,
    fit_polca$llik,
    tolerance = .01
  )

  testthat::expect_equal(
    fit$fit$BIC,
    fit_polca$bic,
    tolerance = .01
  )
})
