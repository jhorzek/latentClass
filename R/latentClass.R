#' @useDynLib latentClass, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

Rcpp::loadModule("NormalModule", TRUE)
Rcpp::loadModule("LCMModule", TRUE)