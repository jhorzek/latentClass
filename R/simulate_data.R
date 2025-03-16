#' simulate_latent_class_data
#'
#' Generates simulated data for a latent class model.
#'
#' @returns data frame with cat_1 and cat_2 (continuous variables) and norm_1 and norm_2 (continuous variables).
#' @importFrom stats rnorm
#' @export
simulate_latent_class_data <- function() {
  # Generate some data for demonstration purposes
  n_samples <- 500
  n_classes <- 3

  # Simulate categorical data
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

  # Simulate normal data
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
    )
  )

  return(cbind(categorical_data, normal_data))
}
