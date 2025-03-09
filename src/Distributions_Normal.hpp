#include <Rcpp.h>
#include <math.h>
#include "Distributions.hpp"

size_t locate_parameter(const std::string& parameter,
                        const std::vector<std::string>& parameter_names) {
  for (size_t i = 0; i < parameter_names.size(); i++) {
    if (parameter.compare(parameter_names.at(i)) == 0) {
      return(i);
    }
  }
  Rcpp::stop("Could not find the parameter: " + parameter + ".");
}

double log_normal(const double x,
                  const double mu,
                  const double sigma){
  
  double ll = -std::log(sigma) -
    0.5 * std::log(2.0*M_PI) -
    0.5 * std::pow((x -  mu) / sigma, 2.0);
  return(ll);
}

double weighted_mean(const std::vector<double>& x,
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

double weighted_n(const std::vector<double>& x,
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

double weighted_standard_deviation(const std::vector<double>& x,
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
  
  struct {
    double mu;
    // We parameterize with log-sigma to prevent negative variances:
    double log_sigma;
  } current_values;
  
  std::vector<bool> free;
  
  void set_pars(const std::vector<std::string>& parameter_names,
                const std::vector<double>& parameter_values){
    if(this->free[0]){
      this->current_values.mu = parameter_values[locate_parameter(this->dist_par_names[0],
                                                                  parameter_names)];
    }
    if(free[1]){
      this->current_values.log_sigma = parameter_values[locate_parameter(this->dist_par_names[1],
                                                                         parameter_names)];
    }
    // Avoid very small variances
    if(std::exp(current_values.log_sigma) < 1e-5){
      current_values.log_sigma = std::log(1e-5);
    }
  }
  
public:
  
  // Constructor
  Normal(const std::vector<double>& data,
         const std::vector<std::string>& dist_par_names,
         std::vector<double> starting_values,
         std::vector<bool> free)
    : Distribution(data, dist_par_names),
      free(free) {
    this->current_values.mu = starting_values[0];
    this->current_values.log_sigma = starting_values[1];
  }
  
  double log_likelihood(const std::vector<std::string>& parameter_names,
                        const std::vector<double>& parameter_values,
                        const std::vector<double>& weights) override {
                          // We assume that the first element in the distribution parameters is the 
                          // mean and the second the log of the standard deviation.
                          
                          this->set_pars(parameter_names,
                                         parameter_values);
                          
                          double n = weighted_n(this->data,
                                                weights);
                          double mean = weighted_mean(this->data,
                                                      weights);
                          double sd = weighted_standard_deviation(this->data,
                                                                  weights,
                                                                  mean);
                          
                          double sigma = std::exp(this->current_values.log_sigma);
                          
                          double ll = n*std::log(2.0*M_PI) + 
                            n*std::log(std::pow(sigma, 2.0)) + 
                            n*std::pow(sd, 2.0)/std::pow(sigma,2) +
                            n * std::pow((mean - this->current_values.mu), 2)/std::pow(sigma,2);
                          ll /= -2.0;
                          return(ll);
                        }
  
  std::vector<double> gradients(const std::vector<std::string>& parameter_names,
                                const std::vector<double>& parameter_values,
                                const std::vector<double>& weights) override {
                                  
                                  this->set_pars(parameter_names,
                                                 parameter_values);
                                  
                                  double n = weighted_n(this->data,
                                                        weights);
                                  double mean = weighted_mean(this->data,
                                                              weights);
                                  double sd = weighted_standard_deviation(this->data,
                                                                          weights,
                                                                          mean);
                                  
                                  // We compute the gradients with respect to 
                                  // each of the parameters in parameter_names.
                                  std::vector<double> gradients(parameter_values.size(), 0.0);
                                  
                                  for(int i = 0; i < parameter_names.size(); i++){
                                    if(parameter_names.at(i).compare(this->dist_par_names.at(0)) == 0){
                                      // derivative with respect to mean
                                      gradients.at(i) = -2.0 * n * (mean - this->current_values.mu)/std::pow(std::exp(this->current_values.log_sigma), 2.0);
                                      gradients.at(i) /= -2.0;
                                    }else if(parameter_names.at(i).compare(this->dist_par_names.at(1)) == 0){
                                      // derivative with respect to standard deviation
                                      gradients.at(i) = 2.0 * n * 
                                        (
                                            1.0 - 
                                              exp(-2.0 * this->current_values.log_sigma) * 
                                              (std::pow(this->current_values.mu, 2.0) -
                                              2.0* this->current_values.mu * mean +
                                              std::pow(sd, 2.0) +
                                              std::pow(mean, 2.0))
                                        );
                                      gradients.at(i) /= -2.0;
                                    }
                                  }
                                  return(gradients);
                                }
  
  std::vector<double> individual_log_likelihood(const std::vector<std::string>& parameter_names,
                                                const std::vector<double>& parameter_values,
                                                const std::vector<double>& weights) override {
                                                  
                                                  if(weights.size() != this->data.size()){
                                                    Rcpp::stop("Weights must be the same size as data.");
                                                  }
                                                  
                                                  this->set_pars(parameter_names,
                                                                 parameter_values);
                                                  
                                                  double sigma = std::exp(this->current_values.log_sigma);
                                                  
                                                  // Compute individual log-likelihood:
                                                  std::vector<double> ind_ll(this->data.size());
                                                  for(int i = 0; i < this->data.size(); i++){
                                                    // For NA values, we basically just return 0. This implies that
                                                    // the current normally distributed variable does not contribute to
                                                    // the individual log-likelihood as the likelihood is
                                                    // multiplied with 1.
                                                    if(isnan(this->data.at(i))){
                                                      ind_ll.at(i) = 0.0;
                                                    }else{
                                                      ind_ll.at(i) = weights.at(i) * log_normal(this->data.at(i), 
                                                                this->current_values.mu, 
                                                                sigma);
                                                    }
                                                  }
                                                  return ind_ll;
                                                }
};