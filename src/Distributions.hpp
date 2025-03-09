#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
#include <Rcpp.h>

// We want to create a vector of distribution classes. To allow for doing
// so, we need a base class that has no templates
class DistributionBase {
public:
  virtual ~DistributionBase() = default;
  
  virtual double log_likelihood(const std::vector<std::string>& parameter_names,
                                const std::vector<double>& parameter_values,
                                const std::vector<double>& weights) = 0;
  
  virtual std::vector<double> individual_log_likelihood(const std::vector<std::string>& parameter_names,
                                                        const std::vector<double>& parameter_values,
                                                        const std::vector<double>& weights) = 0;
  
  virtual std::vector<double> gradients(const std::vector<std::string>& parameter_names,
                                        const std::vector<double>& parameter_values,
                                        const std::vector<double>& weights) = 0;
};

// Based on this class, we now define the templated base class that specific
// distributions inherit from:
template <typename T>
class Distribution : public DistributionBase {
private:
  
public:
  // Member variables
  T data;
  std::vector<std::string> dist_par_names;
  
  // Constructor
  Distribution(const T& data, const std::vector<std::string> dist_par_names)
    : data(data), dist_par_names(dist_par_names) {}
  
  // Virtual destructor
  virtual ~Distribution() = default;
  
  // Pure virtual function for log-likelihood
  virtual double log_likelihood(const std::vector<std::string>& parameter_names,
                                const std::vector<double>& parameter_values,
                                const std::vector<double>& weights) = 0;
  virtual std::vector<double> individual_log_likelihood(const std::vector<std::string>& parameter_names,
                                                        const std::vector<double>& parameter_values,
                                                        const std::vector<double>& weights) = 0;
  virtual std::vector<double> gradients(const std::vector<std::string>& parameter_names,
                                        const std::vector<double>& parameter_values,
                                        const std::vector<double>& weights) = 0;
};

#endif
