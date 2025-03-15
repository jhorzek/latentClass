#ifndef DISTRIBUTIONS_H
#define DISTRIBUTIONS_H
#include <RcppArmadillo.h>

// We want to create a vector of distribution classes. To allow for doing
// so, we need a base class that has no templates
struct model_parameters{
  std::string item_name; 
  arma::mat parameter_values;
  std::vector<std::string> row_names;
  std::vector<std::string> col_names;
};

class DistributionBase {
public:
  virtual ~DistributionBase() = default;
  virtual void maximize_parameters(arma::mat& responsibilities) = 0;
  virtual arma::mat individual_log_likelihood(const std::vector<double>& weights) = 0;
  virtual model_parameters get_parameters() = 0;
};

// Based on this class, we now define the templated base class that specific
// distributions inherit from:
template <typename T>
class Distribution : public DistributionBase {
private:
  
public:
  // Member variables
  T data;
  
  // Constructor
  Distribution(const T& data)
    : data(data) {}
  
  // Virtual destructor
  virtual ~Distribution() = default;
  
  // Pure virtual function for log-likelihood
  virtual void maximize_parameters(arma::mat& responsibilities) = 0;
  virtual arma::mat individual_log_likelihood(const std::vector<double>& weights) = 0;
  virtual model_parameters get_parameters() = 0;
};

#endif
