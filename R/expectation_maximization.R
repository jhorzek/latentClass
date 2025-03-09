
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
                                     conv_crit = 1e-7,
                                     use_ad_gradients = FALSE){
  # initialize responsibilities
  model$update_responsibilities()
  
  pars <- model$get_parameters()
  opt_fun <- function(pars){
    #print(pars)
    model$set_parameters(names(pars),
                         pars)
    return(-2*model$log_likelihood())
  }
  grad_fun <- function(pars){
    model$set_parameters(names(pars),
                         pars)
    return(-2*model$gradients())
  }
  grad_fun(pars)
  ll_old <- model$log_likelihood()
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
    if(use_ad_gradients){
      opt <- optim(par = model$get_parameters(),
                   fn = opt_fun, 
                   gr = grad_fun, 
                   method = "BFGS")
    }else{
      opt <- optim(par = model$get_parameters(),
                   fn = opt_fun,  
                   method = "BFGS")
    }
    model$set_parameters(names(opt$par),
                         opt$par)
    print(opt$par)
    ll_new <- model$log_likelihood()
    if(abs((ll_old - ll_new)/ll_old) < conv_crit){
      converged <- TRUE
      break
    }
    ll_old <- ll_new
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




