#' @useDynLib stantest, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

Rcpp::loadModule("NormalModule", TRUE)