# Export the C++ module to estimate latent classes.
Rcpp::loadModule("LCMModule", TRUE)

#' Estimate a latent class model.
#'
#' latentClass tries to identify subgroups of individuals within the data set. The
#' central assumption is that, once we know the class an individual is part of,
#' all items are independent (local independence). The likelihood of the model
#' is given by:
#'
#' \deqn{L = \prod_{i=1}^N \prod_{c=1}^C \pi_cp(x_{i1}|c)p(x_{i2}|c)\cdots p(x_{iP}|c),}
#'
#' where N is the sample size, C is the number of classes, \eqn{\pi_c} is the
#' probability of class c, \eqn{p(x_{i1}|c)} is the likelihood of observing
#' \eqn{x_{i1}} given that person i is in class c, and P is the number of items.
#'
#' Currently, latentClass supports normal (Gaussian) items and categorical items:
#'
#' \itemize{
#' \item {Gaussian items: The likelihood for normal (Gaussian) items is given by \eqn{\frac{1}{\sqrt{2\pi\sigma^2}}e^{-\frac{(x-\mu)^2}{2\sigma^2}}} (see https://en.wikipedia.org/wiki/Normal_distribution for more details).}
#' \item {Categorical items: The likelihood for categorical items is given by \eqn{p(x) = p_1^{[x=1]}p_2^{[x=2]}\cdots p_k^{[x=k]}} (see https://en.m.wikipedia.org/wiki/Categorical_distribution for more details).}
#' }
#'
#' Parameters are estimated with an Expectation Maximization optimizer. Missing data is
#' handled with full information maximum likelihood (see Vermunt et al., 2013).
#' Additionally, latentClass supports sample weights (see Murphy & Scrucca, 2012).
#'
#' latentClass is one of multiple packages that implement latent class models
#' in R. The implementation of categorical items is heavily inspired by the poLCA
#' package (Drew et al., 2011). Similarly, the implementation of Gaussian items
#' is inspired by mclust (Scrucca, et al., 2023).
#' The implementation of the Expectation Maximization
#' optimizer follows that found in Blume (2002). Overall latentClass is a very
#' vanilla implementation of latent class models, lacking many of the more advanced
#' features of other packages (e.g., no latent class regressions).
#'
#' References:
#'
#' Blume, M. (2002). Expectation maximization: A gentle introduction. Technical University of Munich Institute for Computer Science.
#'
#' Drew A. Linzer, Jeffrey B. Lewis (2011). poLCA: An R Package for Polytomous Variable
#' Latent Class Analysis. Journal of Statistical Software, 42(10), 1-29. URL
#' https://www.jstatsoft.org/v42/i10/
#'
#' Murphy, T. B., & Scrucca, L. (2012). Using Weights in mclust.
#'
#' Scrucca L, Fraley C, Murphy TB, Raftery AE (2023). _Model-Based Clustering,
#' Classification, and Density Estimation Using mclust in R_. Chapman and Hall/CRC. ISBN
#' 978-1032234953, doi:10.1201/9781003277965 <https://doi.org/10.1201/9781003277965>,
#' <https://mclust-org.github.io/book/>
#'
#' Vermunt, J. K., & Magidson, J. (2013). Technical guide for Latent GOLD 5.0:
#' Basic, advanced, and syntax. Belmont, MA: Statistical Innovations Inc.
#'
#' @param data data.frame with the data used for modeling
#' @param categorical object of class Categorical. This object is created with
#' `?latentClass::categorical` and identifies all categorical items in the latent
#' class model
#' @param normal object of class Normal. This object is created with
#' `?latentClass::normal` and identifies all Gaussian items in the latent
#' class model
#' @param n_classes integer defining the number of classes to estimate
#' @param class_probabilities optional: provide starting values for the class
#' probabilities. If provided, a vector of length n_classes that sums to 1.
#' @param sample_weights Not yet implemented: vector with sample weights used in the estimation.
#' @param n_restarts the optimization of latent class models can end up in a local minumum. To try
#' to get to a global minimum, latentClass restarts the the estimation n_restarts times.
#' @param n_cores the number of computer cores to use for the estimation
#' @param opt_settings optimizer settings specified with `latentClass::optimizer_settings`
#'
#' @returns An object of class latentClass with estimates and fit metrics.
#' @importFrom methods is
#' @export
#' @useDynLib latentClass, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom parallel detectCores
#' @importFrom parallel makeCluster
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#'
#' @examples
#' library(latentClass)
#' data <- simulate_latent_class_data()
#' model <- latentClass(data = data,
#'                      # Define the number of classes:
#'                      n_classes = 3,
#'                      # Specify which of the items are categorical:
#'                      categorical = categorical(items = c("cat_1", "cat_2")),
#'                      # Specify which of the items are Gaussian:
#'                      normal = normal(items = c("norm_1", "norm_2")))
#'
#' summary(model)
latentClass <- function(
  data,
  categorical = NULL,
  normal = NULL,
  n_classes,
  class_probabilities = NULL,
  sample_weights = NULL,
  n_restarts = 10,
  n_cores = 2,
  opt_settings = optimizer_settings()
) {
  if (!is.data.frame(data)) stop("data must be a data.frame")
  if (is.null(categorical) & is.null(normal))
    stop("Specify either categorical or normal.")
  if (!is.null(categorical) && !is(categorical, "Categorical"))
    stop(
      "categorical must be of class Categorical (see ?latentClass::categorical)."
    )
  if (!is.null(normal) && !is(normal, "Normal"))
    stop("normal must be of class Normal (see ?latentClass::normal).")
  if (is.null(class_probabilities))
    class_probabilities <- rep(1, n_classes) / n_classes
  if (length(class_probabilities) != n_classes) {
    stop("The length of class_probabilities must be identical to n_classes.")
  }
  if ((sum(class_probabilities) - 1) > 1e-6) {
    stop("class_probabilities must sum to 1.")
  }
  if (!is.null(sample_weights)) {
    if (length(sample_weights) != nrow(data)) {
      stop("sample_weights must have the same length as the data set has rows.")
    }
    if (anyNA(sample_weights)) stop("NAs in sample weights are not allowed")
    if (any(sample_weights < 0))
      warning("Some sample weights are negative. Is this intentional?")
  }

  core_cluster <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(core_cluster)
  on.exit(parallel::stopCluster(core_cluster))

  parallel_results <- foreach::foreach(
    i = seq_len(n_restarts),
    .combine = "c"
  ) %dopar%
    {
      categoricals_initialized <- initialize_categorical(
        data = data,
        categoricals = categorical,
        n_classes = n_classes
      )

      normals_initialized <- initialize_normal(
        data = data,
        normals = normal,
        n_classes = n_classes
      )

      # set up the latent class model
      model <- LCMR$new(class_probabilities, nrow(data))
      if (!is.null(sample_weights)) {
        model$set_sample_weights(sample_weights)
      }

      for (i in seq_len(length(categoricals_initialized$items))) {
        model$add_categorical(
          categoricals_initialized$items[[i]],
          factor_to_index(data[[categoricals_initialized$items[[
            i
          ]]]]),
          categoricals_initialized$starting_values[[i]]
        )
      }
      for (i in seq_len(length(normals_initialized$items))) {
        model$add_normal(
          normals_initialized$items[[i]],
          data[[normals_initialized$items[[i]]]],
          # means
          normals_initialized$starting_values[[i]][1, ],
          # standard deviations
          normals_initialized$starting_values[[i]][2, ],
          normals_initialized$sd_equal[i]
        )
      }

      # Optimize
      timing <- system.time({
        converged <- model$expectation_maximization(
          opt_settings$max_iter,
          opt_settings$convergence_criterion
        )
      })
      result <- finalize_estimates(
        data = data,
        sample_weights = sample_weights,
        model = model,
        categoricals_initialized = categoricals_initialized,
        normals_initialized = normals_initialized,
        converged = converged,
        timing
      )
      class(result) <- "latentClass"
      ret <- list(result)
      names(ret) <- paste0("iteration_", i)
      ret
    }
  best_fit <- which.max(sapply(
    parallel_results,
    function(x) x$fit$`log-Likelihood`
  ))[1]

  return(parallel_results[[best_fit]])
}

#' print.latentClass
#'
#' Print the results of a latent class model.
#'
#' @param x A latent class model object.
#' @param ... not used
#' @returns nothing
#' @export
print.latentClass <- function(x, ...) {
  cat("#### Latent Class Model Results #####\n\n")

  # Describe model settings
  cat("Model settings:\n")
  cat(paste0(rep("-", nchar("Model settings")), collapse = ""), "\n")
  max_char <- max(
    nchar("- Number of classes: "),
    nchar("- Categorical variables: "),
    nchar("- Normal variables: ")
  )
  cat(paste0(
    "- Number of classes: ",
    paste0(rep(" ", max_char - nchar("- Number of classes: ")), collapse = ""),
    x$n_classes,
    "\n"
  ))
  cat(paste0(
    "- Categorical variables: ",
    paste0(
      rep(
        " ",
        max_char - nchar("- Categorical variables: ")
      ),
      collapse = ""
    ),
    paste0(x$categorical_items, collapse = ", "),
    "\n"
  ))
  cat(paste0(
    "- Normal variables: ",
    paste0(rep(" ", max_char - nchar("- Normal variables: ")), collapse = ""),
    paste0(x$normal_items, collapse = ", "),
    "\n\n"
  ))

  # Model estimation
  cat("Estimation:\n")
  cat(paste0(rep("-", nchar("Estimation")), collapse = ""), "\n")
  cat("- Model converged:", ifelse(x$converged, "Yes", "No"), "\n")
  cat("- Estimation time:", x$processing_time["elapsed"], "seconds", "\n")
  cat("\n")

  # Fit measures
  cat("Fit measures:\n")
  cat(paste0(rep("-", nchar("Fit measures")), collapse = ""), "\n")
  fit_nchar <- max(sapply(names(x$fit), nchar))
  for (fm in names(x$fit)) {
    cat(paste0(
      "- ",
      fm,
      ": ",
      paste0(rep(" ", fit_nchar - nchar(fm)), collapse = ""),
      ifelse(x$fit[[fm]] > 0, " ", ""),
      round(x$fit[[fm]], 3),
      "\n"
    ))
  }
  cat("\n")

  cat("Estimates:\n")
  cat(paste0(rep("-", nchar("Estimates")), collapse = ""), "\n")
  for (item in names(x$estimates)) {
    cat(paste0(item, ":\n"))
    print(x$estimates[[item]])
    cat("\n")
  }
}

#' summary.latentClass
#'
#' Print the results of a latent class model.
#'
#' @param object A latent class model object.
#' @param ... not used
#' @returns nothing
#' @export
summary.latentClass <- function(object, ...) {
  print(object)
}

#' Specify the items that follow a normal (Gaussian) distribution.
#'
#' For all items specified here, latentClass will use a normal distribution
#' density as likelihood.
#'
#' @param items character vector with the names of the items
#' @param sd_equal if set to TRUE, the standard deviations for each item
#' will be set to equal across classes. Alternatively, a vector can be passed
#' that specifies for each of the items if the standard deviations should be equal
#' across classes
#' @param starting_values if set to "random", starting values are generated randomly.
#' Alternatively, a list with matrices holding starting values for each item can be provided. Each of the
#' starting values matrices must have 2 rows (starting values for mean and standard deviation)
#' and as many columns as there are classes
#'
#' @return An object of class Normal.
#' @export
#'
#' @examples
#' library(latentClass)
#' normal(items = c("x1", "x2"))
#' # Allow for class-specific standard deviations:
#' normal(items = c("x1", "x2"),
#'        sd_equal = FALSE)
#' # Allow for class-specific standard deviations for x1, but not x2:
#' normal(items = c("x1", "x2"),
#'        sd_equal = c(TRUE, FALSE))
#' # Provide starting values (assuming there are 3 classes):
#' normal(items = c("x1", "x2"),
#'        starting_values =
#'        list("x1" = matrix(c(-3, 0, 3,
#'                              1, 1, 1),
#'                          nrow = 2,
#'                          byrow = TRUE,
#'                          dimnames = list(c("mean", "sd"),
#'                                          c("class_1", "class_2", "class_3"))),
#'             "x2" = matrix(c(-1, 0, 1,
#'                              1, 1, 1),
#'                          nrow = 2,
#'                          byrow = TRUE,
#'                          dimnames = list(c("mean", "sd"),
#'                                          c("class_1", "class_2", "class_3")))))
normal <- function(items, sd_equal = TRUE, starting_values = "random") {
  ret <- list(
    items = items,
    sd_equal = sd_equal,
    starting_values = starting_values
  )
  class(ret) <- "Normal"
  return(ret)
}

#' @importFrom stats quantile
#' @importFrom stats runif
#' @importFrom stats sd
#' @noRd
initialize_normal <- function(data, normals, n_classes) {
  if (is.null(normals)) {
    return(list(
      items = c()
    ))
  }

  check_distribution(data, normals)

  if (length(normals$sd_equal) == 1) {
    normals$sd_equal <- rep(normals$sd_equal, length(normals$items))
    names(normals$sd_equal) <- normals$items
  } else {
    if (length(normals$sd_equal) != length(normals$items))
      stop("Please specify one sd_equal setting per item.")
    names(normals$sd_equal) <- normals$items
  }

  if (is.list(normals$starting_values)) {
    if (
      length(union(names(normals$starting_values), normals$items)) !=
        length(normals$items)
    )
      stop("Please provide starting values for all items.")
    for (starting_values in normals$starting_values) {
      if (!is.matrix(starting_values))
        stop("starting_values must be a list of matrices.")
      if (!is.numeric(starting_values))
        stop("starting_values must be a numeric matrix.")
      if (nrow(starting_values) != 2)
        stop(
          "starting_values of a normal distribution must have two rows: The first row ",
          "has the starting values for the means and the second the starting values for ",
          "the standard deviations."
        )
      if (ncol(starting_values) != n_classes)
        stop(
          "starting_values of a normal distribution must as many columns as there are classes."
        )
    }
  } else if (normals$starting_values == "random") {
    starting_values <- list()
    for (item in normals$items) {
      starting_values[[item]] <- matrix(0, nrow = 2, ncol = n_classes)
      # We want to make sure that we are not fully off here, so we create starting
      # values randomly based on the distribution
      starting_values[[item]][1, ] <- stats::runif(
        n = n_classes,
        min = stats::quantile(data[[item]], .05, na.rm = TRUE),
        max = stats::quantile(data[[item]], .95, na.rm = TRUE)
      )
      if (normals$sd_equal[item]) {
        starting_values[[item]][2, ] <- rep(
          .5 * stats::sd(data[[item]], na.rm = TRUE),
          n_classes
        )
      } else {
        starting_values[[item]][2, ] <- stats::runif(
          n = n_classes,
          min = .1 * stats::sd(data[[item]], na.rm = TRUE),
          max = 4 * stats::sd(data[[item]], na.rm = TRUE)
        )
      }
    }
    normals$starting_values <- starting_values
  } else {
    stop("Unknown setting for the starting values.")
  }
  return(normals)
}

#' Specify the items that follow a normal categorical distribution.
#'
#' @param items character vector with the names of the items
#' @param starting_values if set to "random", starting values are generated randomly.
#' Alternatively, a list with matrices holding starting values for each item can be provided. Each of the
#' starting values matrices must have as many rows as the item has levels. Additionally,
#' the columns of each matrix must sum to 1.
#'
#' @return An object of class Categorical.
#' @export
#'
#' @examples
#' library(latentClass)
#' categorical(items = c("x1", "x2"))
#' # Provide starting values (assuming there are 3 classes and x1 has 2 levels,
#' # while c2 has 3 levels):
#' categorical(items = c("x1", "x2"),
#'        starting_values =
#'        list("x1" = matrix(c(.3, .5, .8,
#'                             .7, .5, .2),
#'                          nrow = 2,
#'                          byrow = TRUE,
#'                          dimnames = list(c("level_1", "level_2"),
#'                                          c("class_1", "class_2", "class_3"))),
#'             "x2" = matrix(c(.1, .7, .2,
#'                             .8, .2, .5,
#'                             .1, .1, .3),
#'                          nrow = 3,
#'                          byrow = TRUE,
#'                          dimnames = list(c("level_1", "level_2", "level_3"),
#'                                          c("class_1", "class_2", "class_3"))))
#' )
categorical <- function(items, starting_values = "random") {
  ret <- list(
    items = items,
    starting_values = starting_values
  )
  class(ret) <- "Categorical"
  return(ret)
}

#' initialize_categorical
#'
#' Initializes the setup for categorical items (e.g., set up
#' starting values).
#'
#' @param data data set
#' @param categoricals categorical setup provided by the user
#' @param n_classes number of classes
#' @returns list with setup
#' @importFrom stats quantile
#' @importFrom stats runif
#' @noRd
initialize_categorical <- function(data, categoricals, n_classes) {
  if (is.null(categoricals)) {
    return(list(
      items = c()
    ))
  }

  check_distribution(data, categoricals)

  if (is.list(categoricals$starting_values)) {
    if (
      length(union(names(categoricals$starting_values), categoricals$items)) !=
        length(categoricals$items)
    )
      stop("Please provide starting values for all items.")
    for (starting_values in categoricals$starting_values) {
      if (!is.matrix(starting_values))
        stop("starting_values must be a list of matrices.")
      if (!is.numeric(starting_values))
        stop("starting_values must be a numeric matrix.")

      if (any(abs(colSums(starting_values) - 1) < .0001))
        stop("Columns in starting_values must sum to 1.")

      if (!is.factor(data[[categoricals$items]]))
        stop(paste0(item, " must be a factor."))
      item_levels <- levels(data[[categoricals$items]])

      if (nrow(starting_values) != length(item_levels))
        stop(
          "starting_values of a categorical item must have as many rows as there are levels in the factor."
        )
      if (ncol(starting_values) != n_classes)
        stop(
          "starting_values of a normal distribution must as many columns as there are classes."
        )
    }
  } else if (categoricals$starting_values == "random") {
    starting_values <- list()
    for (item in categoricals$items) {
      if (!is.factor(data[[item]]))
        stop("Item ", item, " must be a factor to be used as categorical.")
      item_levels <- levels(data[[item]])
      starting_values[[item]] <- matrix(
        stats::runif(n = length(item_levels) * n_classes, min = .1, 1),
        nrow = length(item_levels),
        ncol = n_classes
      )
      # Normalize to sum to 1:
      for (cl in 1:n_classes) {
        starting_values[[item]][, cl] <- starting_values[[item]][, cl] /
          sum(starting_values[[item]][, cl])
      }
    }
    categoricals$starting_values <- starting_values
  } else {
    stop("Unknown setting for the starting values.")
  }
  return(categoricals)
}

#' check_distribution
#'
#' Check distribution ensures that the user passed the correct settings
#' in the distribution functions.
#' @param data data set
#' @param dist the distribtion object that should be checked
#' @returns nothing
#' @noRd
check_distribution <- function(data, dist) {
  if (length(dist$items) == 0) invisible(return())

  if (!is.character(dist$items))
    stop("The items of a distribution must be a character vector.")

  if (length(setdiff(dist$items, colnames(data))) > 0) {
    stop(paste0(
      "Could not find ",
      paste0(setdiff(dist$items, colnames(data)), collapse = ", "),
      " in the data."
    ))
  }
}

#' factor_to_index
#'
#' For categorical variables, users may have supplied numerical or string
#' variables. We want to change this to integers starting from zero to
#' 1 minus the number of levels of the factor. This will allow us to use these
#' integers are indices in the computations in C++.
#' @param item the categorical item
#' @returns vector with indices
#' @noRd
factor_to_index <- function(item) {
  item_levels <- levels(item)
  # C++ has no NA integer value -> we use -99 for NA
  index <- rep(-99, length(item))
  index[!is.na(item)] <- sapply(
    item[!is.na(item)],
    # for C++, we need indices starting at 0:
    function(x) which(item_levels == x) - 1
  )
  # for C++, we need indices starting at 0:
  return(index)
}

#' finalize_estimates
#'
#' Combines model parameters and model information to be returned to the user
#' @param data data set
#' @param sample_weights sample weights
#' @param model estimated model
#' @param categoricals_initialized list with info on categorical items
#' @param normals_initialized list with info on normal items
#' @param converged boolean indiccating if the model converged
#' @param timing info on the time it took to estimate the model
#' @returns list with model info
#' @noRd
finalize_estimates <- function(
  data,
  sample_weights,
  model,
  categoricals_initialized,
  normals_initialized,
  converged,
  timing
) {
  ll <- model$log_likelihood()
  class_probabilities <- model$get_class_probabilities()
  names(class_probabilities) <- paste0("class_", 1:length(class_probabilities))
  pars <- model$get_parameters()
  # for categorical items, we have to replace the internal factor levels in the
  # rownames with those from the items
  for (item in categoricals_initialized$items) {
    lvls <- levels(data[[item]])
    names(lvls) <- as.character((1:length(lvls)) - 1)
    rownames(pars[[item]]) <- lvls[rownames(pars[[item]])]
  }

  # Count the number of parameters
  n_parameters <- length(class_probabilities) - 1 # if we know c-1 class probabilities, we also
  # now the last one.
  for (item in categoricals_initialized$items) {
    lvls <- levels(data[[item]])
    # if we know l-1 level probabilities, we also
    # now the last one.
    n_parameters <- n_parameters + ((length(lvls) - 1) * model$get_n_classes())
  }
  for (item in normals_initialized$items) {
    n_parameters <- n_parameters +
      (
        # means
        model$get_n_classes() +
          # standard deviations
          model$get_n_classes()^(!normals_initialized$sd_equal[item])
      )
  }

  fit <- fit_measures(
    ll = ll,
    n = nrow(data),
    n_parameters = n_parameters
  )

  res <- list(
    n_classes = model$get_n_classes(),
    class_probabilities = class_probabilities,
    categorical_items = categoricals_initialized$items,
    normal_items = normals_initialized$items,
    estimates = pars,
    fit = fit,
    sample_weights,
    converged = converged,
    processing_time = timing
  )

  return(res)
}

#' fit_measures
#'
#' Computes fit measures for the estimated model
#' @param ll log-likelihood
#' @param n number of persons
#' @param n_parameters number of parameters
#' @returns list with fit measures
#' @noRd
fit_measures <- function(ll, n, n_parameters) {
  # For the computation of the fit measures we follow the implementaiton in mclust, where
  # the sample size instead of the weighted n is used:
  # https://github.com/cran/mclust/blob/65e2a1c0538807f5e52ade030f12f5af4c1bc746/R/weights.R#L67C19-L67C33
  bic <- -2 * ll + log(n) * n_parameters
  aic <- -2 * ll + 2 * n_parameters

  return(list(
    "Prameters" = unname(n_parameters),
    "Observations" = unname(n),
    "log-Likelihood" = unname(ll),
    BIC = unname(bic),
    AIC = unname(aic)
  ))
}

#' Fine tune the optimizer settings.
#'
#' @param max_iter maximum number of iterations
#' @param convergence_criterion cut-off value. If the relative change in the
#' log-likelihood from one iteration to the next falls below this value, the optimization
#' stops and the model is seen as converged.
#'
#' @return list with optimizer settings
#' @export
#'
#' @examples
#' library(latentClass)
#' optimizer_settings()
optimizer_settings <- function(max_iter = 1000, convergence_criterion = 1e-10) {
  return(
    list(max_iter = max_iter, convergence_criterion = convergence_criterion)
  )
}
