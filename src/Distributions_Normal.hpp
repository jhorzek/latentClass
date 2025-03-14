#ifndef NORMAL_H
#define NORMAL_H
#include <math.h>
#include "Distributions.hpp"

inline double log_normal(const double x,
                         const double mu,
                         const double sigma){
  
  double ll = -std::log(sigma) -
    0.5 * std::log(2.0*M_PI) -
    0.5 * std::pow((x -  mu) / sigma, 2.0);
  return(ll);
}

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

// Derived class for normal distribution
class Normal : public Distribution<std::vector<double>> {
private:
  
  std::string item_name;
  // parameter struct holding the means and standard deviations for 
  // all classes
  struct {
    std::vector<double> means;
    std::vector<double> sds;
    bool sd_equal;
  } parameters;
  
  int n_classes;
  
public:
  
  // Constructor
  Normal(std::string item_name,
         const std::vector<double>& data,
         std::vector<double> means,
         std::vector<double> sds,
         bool sd_equal)
    : Distribution(data) {
    if(means.size() != sds.size())
      Rcpp::stop("Means and sds must be of the same length");
    this->item_name = item_name;
    this->parameters.means = means;
    this->parameters.sds = sds;
    this->parameters.sd_equal = sd_equal;
    this->n_classes = means.size();
  }
  
  void maximize_parameters(arma::mat& responsibilities) override {
    // We have to compute the means and standard deviations for each of the 
    // classes.
    
    double squared_residual_sum = 0.0;
    double n = 0.0;
    for(int cl = 0; cl < this->n_classes; cl++){
      if(!this->parameters.sd_equal){
        squared_residual_sum = 0.0;
        n = 0.0;
      }
      
      // Convert the Armadillo vector to a std::vector<double>
      std::vector<double> current_responsibilities(responsibilities.col(cl).begin(), 
                                                   responsibilities.col(cl).end());
      
      this->parameters.means.at(cl) = weighted_mean(this->data,
                                current_responsibilities);
      
      for(int i = 0; i < this->data.size(); i++){
        squared_residual_sum += current_responsibilities.at(i)*
          (std::pow(this->data.at(i) - 
          this->parameters.means.at(cl), 2.0));
        n += current_responsibilities.at(i);
      }
      if(!this->parameters.sd_equal){
        this->parameters.sds.at(cl) = std::sqrt(squared_residual_sum / n);
      }
    }
    
    if(this->parameters.sd_equal){
      // If the standard deviations are equal, we simply compute the sd over all classes.
      // The main thing here is that each data point is compared to the class specific
      // mean when forming the squared residuals, which is what we have done above.
      for(int cl = 0; cl < this->n_classes; cl++){
        this->parameters.sds.at(cl) = std::sqrt(squared_residual_sum / n);
      }
    }
  }
  
  arma::mat individual_log_likelihood(const std::vector<double>& weights) override {
    
    if(weights.size() != this->data.size()){
      Rcpp::stop("Weights must be of the same size as the data.");
    }
    // Returns the log-likelihood for every person (rows) and class (columns).
    // Note that the log-likelihood is unweighted.
    arma::mat log_lik(this->data.size(), this->n_classes);
    
    for(std::size_t i = 0; i < this->data.size(); i++){
      for(std::size_t cl = 0; cl < this->n_classes; cl++){
        log_lik(i, cl) = weights.at(i) * 
          log_normal(this->data.at(i),
                     this->parameters.means.at(cl),
                     this->parameters.sds.at(cl));
      }
    }
    
    return log_lik;
  }
  
  model_parameters get_parameters() override {
    
    model_parameters params;
    params.parameter_values.reshape(2, this->n_classes);
    params.item_name = this->item_name;
    params.row_names = {"mean", "sd"};
    
    for(int cl = 0; cl < this->n_classes; cl++){
      params.col_names.push_back("class_" + std::to_string(cl+1));
      
      params.parameter_values(0,cl) = this->parameters.means.at(cl);
      params.parameter_values(1,cl) = this->parameters.sds.at(cl);
    }
    
    return(params);
  } 
  
};

#endif