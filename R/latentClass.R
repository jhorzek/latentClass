#' @useDynLib latentClass, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

Rcpp::loadModule("LCMModule", TRUE)


latentClass <- function(data,
                        categorical = NULL,
                        normal = NULL,
                        n_classes,
                        class_probabilities = NULL,
                        opt_settings = optimizer_settings()){
  
  if(!is.data.frame(data))
    stop("data must be a data.frame")
  if(is.null(categorical) & is.null(normal))
    stop("Specify either categorical or normal.")
  if(!is.null(categorical) && !is(categorical, "Categorical"))
    stop("categorical must be of class Categorical (see ?latentClass::categorical).")
  if(!is.null(normal) && !is(normal, "Normal"))
    stop("normal must be of class Normal (see ?latentClass::normal).")
  if(is.null(class_probabilities))
    class_probabilities <- rep(1, n_classes)/n_classes
  if(length(class_probabilities) != n_classes){
    stop("The length of class_probabilities must be identical to n_classes.")
  }
  if((sum(class_probabilities) - 1) > 1e-6){
    stop("class_probabilities must sum to 1.")
  }
  
  categoricals_initialized <- initialize_categorical(data = data, 
                                                     categoricals = categorical, 
                                                     n_classes = n_classes)
  
  normals_initialized <- initialize_normal(data = data, 
                                           normals = normal, 
                                           n_classes = n_classes)
  
  # set up the latent class model
  model <- latentClass:::LCMR$new(class_probabilities,
                                  nrow(data))
  for(i in seq_len(length(categoricals_initialized$items))){
    model$add_categorical(categoricals_initialized$items[[i]],
                          factor_to_index(data[[categoricals_initialized$items[[i]]]]),
                          categoricals_initialized$starting_values[[i]])
  }
  for(i in seq_len(length(normals_initialized$items))){
    model$add_normal(normals_initialized$items[[i]],
                     data[[normals_initialized$items[[i]]]],
                     # means
                     normals_initialized$starting_values[[i]][1,],
                     # standard deviations
                     normals_initialized$starting_values[[i]][2,],
                     normals_initialized$sd_equal[i])
  }
  
  # Optimize
  converged <- model$expectation_maximization(opt_settings$max_iter,
                                              opt_settings$convergence_criterion)
  result <- finalize_estimates(data = data, 
                               model = model, 
                               categoricals_initialized = categoricals_initialized, 
                               normals_initialized = normals_initialized, 
                               converged = converged)
  class(result) <- "latentClass"
  return(result)
}

normal <- function(items,
                   sd_equal,
                   starting_values = "random"){
  ret <- list(
    items = items,
    sd_equal = sd_equal,
    starting_values = starting_values
  )
  class(ret) <- "Normal"
  return(ret)
}

initialize_normal <- function(data,
                              normals,
                              n_classes){
  if(is.null(normals)){
    return(list(
      items = c()
    ))
  }
  
  check_distribution(data,
                     normals)
  
  if(length(normals$sd_equal) == 1){
    normals$sd_equal <- rep(normals$sd_equal, length(normals$items))
    names(normals$sd_equal) <- items
  }else{
    if(length(normals$sd_equal) != length(normals$items))
      stop("Please specify one sd_equal setting per item.")
    names(normals$sd_equal) <- items
  }
  
  
  if(is.list(normals$starting_values)){
    if(length(union(names(normals$starting_values), normals$items)) != length(normals$items))
      stop("Please provide starting values for all items.")
    for(starting_values in normals$starting_values){
      if(!is.matrix(starting_values))
        stop("starting_values must be a list of matrices.")
      if(!is.numeric(starting_values))
        stop("starting_values must be a numeric matrix.")
      if(nrow(starting_values) != 2)
        stop("starting_values of a normal distribution must have two rows: The first row ",
             "has the starting values for the means and the second the starting values for ",
             "the standard deviations.")
      if(ncol(starting_values) != n_classes)
        stop("starting_values of a normal distribution must as many columns as there are classes.")
    }
  }else if(normals$starting_values == "random"){
    starting_values <- list()
    for(item in normals$items){
      starting_values[[item]] <- matrix(0, 
                                        nrow = 2, 
                                        ncol = n_classes)
      # We want to make sure that we are not fully off here, so we create starting
      # values randomly based on the distribution
      starting_values[[item]][1, ] <- runif(n = n_classes,
                                            min = quantile(data[[item]], .05, na.rm = TRUE),
                                            max = quantile(data[[item]], .95, na.rm = TRUE))
      if(normals$sd_equal[item]){
        starting_values[[item]][2, ] <- rep(.5*sd(data[[item]], na.rm = TRUE), n_classes)
      }else{
        starting_values[[item]][2, ] <- runif(n = n_classes,
                                              min = .1*sd(data[[item]], na.rm = TRUE),
                                              max = 4*sd(data[[item]], na.rm = TRUE))
      }
    }
    normals$starting_values <- starting_values
  }else{
    stop("Unknown setting for the starting values.")
  }
  return(normals)
}

categorical <- function(items,
                        starting_values = "random"){
  ret <- list(
    items = items,
    starting_values = starting_values
  )
  class(ret) <- "Categorical"
  return(ret)
}

initialize_categorical <- function(data,
                                   categoricals,
                                   n_classes){
  if(is.null(categoricals)){
    return(list(
      items = c()
    ))
  }
  
  check_distribution(data,
                     categoricals)
  
  if(is.list(categoricals$starting_values)){
    if(length(union(names(categoricals$starting_values), categoricals$items)) != length(categoricals$items))
      stop("Please provide starting values for all items.")
    for(starting_values in categoricals$starting_values){
      if(!is.matrix(starting_values))
        stop("starting_values must be a list of matrices.")
      if(!is.numeric(starting_values))
        stop("starting_values must be a numeric matrix.")
      
      if(any(abs(colSums(starting_values) - 1) < .0001))
        stop("Columns in starting_values must sum to 1.")
      
      if(!is.factor(data[[categoricals$items]]))
        stop(paste0(item, " must be a factor."))
      item_levels <- levels(data[[categoricals$items]])
      
      if(nrow(starting_values) != length(item_levels))
        stop("starting_values of a categorical item must have as many rows as there are levels in the factor.")
      if(ncol(starting_values) != n_classes)
        stop("starting_values of a normal distribution must as many columns as there are classes.")
    }
  }else if(categoricals$starting_values == "random"){
    starting_values <- list()
    for(item in categoricals$items){
      if(!is.factor(data[[item]]))
        stop("Item ", item, " must be a factor to be used as categorical.")
      item_levels <- levels(data[[item]])
      starting_values[[item]] <- matrix(runif(n = length(item_levels)*n_classes, min = .1, 1), 
                                        nrow = length(item_levels), 
                                        ncol = n_classes)
      # Normalize to sum to 1:
      for(cl in 1:n_classes){
        starting_values[[item]][,cl] <- starting_values[[item]][,cl]/sum(starting_values[[item]][,cl])
      }
    }
    categoricals$starting_values <- starting_values
  }else{
    stop("Unknown setting for the starting values.")
  }
  return(categoricals)
}

check_distribution <- function(data,
                               dist){
  if(length(dist$items) == 0)
    invisible(return())
  
  if(!is.character(dist$items))
    stop("The items of a distribution must be a character vector.")
  
  if(length(setdiff(dist$items, colnames(data))) > 0){
    stop(paste0("Could not find ", paste0(setdiff(dist$items, colnames(data)), collapse = ", "),
                " in the data."))
  }
  
}

factor_to_index <- function(item){
  
  index <- sapply(item, function(x) which(levels(item) == x))
  # for C++, we need indices starting at 0:
  return(index-1)
}

finalize_estimates <- function(data,
                               model,
                               categoricals_initialized,
                               normals_initialized,
                               converged){
  ll <- model$log_likelihood()
  class_probabilities <- model$get_class_probabilities()
  names(class_probabilities) <- paste0("class_", 1:length(class_probabilities))
  pars <- model$get_parameters()
  # for categorical items, we have to replace the internal factor levels in the 
  # rownames with those from the items
  for(item in categoricals_initialized$items){
    lvls <- levels(data[[item]])
    names(lvls) <- as.character((1:length(lvls)) - 1)
    rownames(pars[[item]]) <- lvls[rownames(pars[[item]])]
  }
  
  # Count the number of parameters
  n_parameters <- 0
  for(item in categoricals_initialized$items){
    lvls <- levels(data[[item]])
    n_parameters <- n_parameters + (length(lvls) * model$get_n_classes())
  }
  for(item in normals_initialized$items){
    n_parameters <- n_parameters  + n_parameters^normals_initialized$sd_equal[item]
  }
  
  fit <- fit_measures(ll = ll,
                      n = nrow(data),
                      n_parameters = n_parameters)
  
  res <- list(
    n_classes = model$get_n_classes(),
    class_probabilities = class_probabilities,
    estimates = pars,
    fit = fit,
    converged = converged
  )
  
  return(res)
}

fit_measures <- function(ll,
                         n,
                         n_parameters){
  bic <- -2*ll + log(n)*n_parameters
  aic <- -2*ll + 2*n_parameters
  
  return(list(
    log_likelihood = ll,
    BIC = bic,
    AIC = aic
  ))
}

optimizer_settings <- function(max_iter = 1000,
                               convergence_criterion = 1e-10){
  return(
    list(max_iter = max_iter,
         convergence_criterion = convergence_criterion)
  )
}