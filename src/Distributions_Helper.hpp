#ifndef DISTHELPER_H
#define DISTHELPER_H
#include <math.h>

inline double weighted_mean(const std::vector<double>& x,
                            const std::vector<double>& weights){
  if(x.size() != weights.size()){
    Rcpp::stop("x must have the same size as weights.");
  }
  double mean = 0.0;
  double n = 0.0;
  for(int i = 0; i < x.size(); i++){
    if(std::isnan(x[i]))
      continue;
    mean += weights[i]*x[i];
    n += weights[i];
  }
  mean /= n;
  return(mean);
}

inline double weighted_n(const std::vector<double>& x,
                         const std::vector<double>& weights){
  if(x.size() != weights.size()){
    Rcpp::stop("x must have the same size as weights.");
  }
  double n = 0.0;
  for(int i = 0; i < x.size(); i++){
    if(std::isnan(x[i]))
      continue;
    n += weights[i];
  }
  return(n);
}

inline double weighted_standard_deviation(const std::vector<double>& x,
                                          const std::vector<double>& weights, 
                                          const double w_mean){
  if(x.size() != weights.size()){
    Rcpp::stop("x must have the same size as weights.");
  }
  double sd = 0.0;
  double n = 0.0;
  for(int i = 0; i < x.size(); i++){
    if(std::isnan(x[i]))
      continue;
    sd += weights[i]*pow((x[i] - w_mean), 2.0);
    n += weights[i];
  }
  sd /= n;
  sd = std::sqrt(sd);
  return(sd);
}

#endif