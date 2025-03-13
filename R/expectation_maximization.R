
#' @examples
#' x <- c(rnorm(5, -3, 1),
#' rnorm(5, 3, 1))
#' model <- stantest:::LCMR$new(c(.5,.5),
#'                              length(x))
#' model$get_n_classes()
#' model$get_class_probabilities()
#' model$get_n_persons()
#' model$add_normal(1,
#'                  x,
#'                  c("mu_1", "log_sd_1"),
#'                  c(-3,0))
#' model$add_normal(2,
#'                  x,
#'                  c("mu_2", "log_sd_1"),
#'                  c(3,0))
#' model$get_responsibilities()
#' 
#' expectation_maximization(model)
#' round(model$get_responsibilities(), 3)
#' model$get_class_probabilities()
#' model$get_parameters()
expectation_maximization <- function(model,
                                     max_iter = 100,
                                     conv_crit = 1e-5){
  # initialize responsibilities
  model$update_responsibilities()
  
  pars_old <- model$get_parameters()
  # ll_old <- model$expected_log_likelihood(names(pars),
  #                                         pars)
  n_iter <- 0
  
  pb <- utils::txtProgressBar(min = 0, max = max_iter, style = 3)
  for(i in 1:max_iter){
    n_iter <- n_iter + 1
    #### Expectation ####
    model$update_responsibilities()
    # responsibilities <- model$get_responsibilities()
    
    if(i > 1) # in the first iteration, we are using the starting values for the class probabilities
      model$update_class_probabilities()
    
    #### Maximization ####
    model$maximize()
    pars_new <- model$get_parameters()
    
    # ll_new <- model$expected_log_likelihood(names(opt$par),
    #                                         opt$par)
    if(max(abs((pars_old - pars_new)/pars_old)) < conv_crit){
      converged <- TRUE
      break
    }
    pars_old <- pars_new
    if(i == max_iter){
      converged <- FALSE
      warning("EM algorithm did not converge.")
    }
    utils::setTxtProgressBar(pb = pb, value = i)
  }
  
  # final expectations update
  model$update_responsibilities()
  model$update_class_probabilities()
  
  return(list(model = model,
              converged = converged,
              n_iterations = n_iter))
}




