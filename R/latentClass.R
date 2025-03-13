#' @useDynLib latentClass, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

Rcpp::loadModule("LCMModule", TRUE)