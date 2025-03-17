#ifndef MULTINOMIAL_H
#define MULTINOMIAL_H
#include <math.h>
#include "Distributions.hpp"
#include "Distributions_Helper.hpp"

// Derived class for normal distribution
class Categorical : public Distribution<std::vector<int>> {
private:
  
  std::string item_name;
  // parameters holding the probabilities for 
  // all classes
  arma::mat parameters;
  
  int n_classes;
  int n_levels; // levels of the factor
  
public:
  
  // Constructor
  Categorical(std::string item_name,
              // the data is in integers. Each integer indicates which option
              // of the factor was selected starting with 0 and ending with options-1.
              // This way, we can use the integers directly as indices for the parameters.
              const std::vector<int>& data,
              arma::mat parameters)
    : Distribution(data) {
    
    this->item_name = item_name;
    this->n_classes = parameters.n_cols;
    this->n_levels = parameters.n_rows;
    this->parameters = parameters;
    for(int i = 0; i < data.size(); i++){
      if((data.at(i) < 0) || (data.at(i) > this->n_levels-1)){
        Rcpp::stop("The data should be indices starting with zero and ending with n_levels - 1");
      }
    }
  }
  
  void maximize_parameters(arma::mat& responsibilities, std::vector<double> sample_weights) override {
    // We have to compute the probabilities for each of the options in each of the
    // classes.
    if(sample_weights.size() != this->data.size()){
      Rcpp::stop("sample_weights must be of the same size as the data.");
    }

    std::vector<double> class_n(this->n_classes, 0.0);

    this->parameters.fill(0.0);
    
    for(int cl = 0; cl < this->n_classes; cl++){
      for(int i = 0; i < this->data.size(); i ++){
        // in case of nan: does not contribute
        if(std::isnan(this->data.at(i)))
        continue;
        // because we implemented the data as indices, we can simply access 
        // the class with the data as index.
        this->parameters(this->data.at(i), cl) += sample_weights.at(i) * responsibilities(i,cl);
        class_n.at(cl) += sample_weights.at(i) * responsibilities(i,cl);
      }
    }
    // Now we just have to divide each column of parameters by the class-specific
    // sample size.
    for(int cl = 0; cl < this->n_classes; cl++){
      this->parameters.col(cl) /= class_n.at(cl);
    }
  }
  
  arma::mat individual_log_likelihood(const std::vector<double>& sample_weights) override {
    
    if(sample_weights.size() != this->data.size()){
      Rcpp::stop("sample_weights must be of the same size as the data.");
    }
    // Returns the log-likelihood for every person (rows) and class (columns).
    arma::mat log_lik(this->data.size(), this->n_classes);
    
    for(std::size_t i = 0; i < this->data.size(); i++){
      for(std::size_t cl = 0; cl < this->n_classes; cl++){
        // in case of nan: does not contribute
        if(std::isnan(this->data.at(i)))
        continue;
        log_lik(i, cl) = sample_weights.at(i) * std::log(this->parameters(this->data.at(i), cl));
      }
    }
    
    return log_lik;
  }
  
  model_parameters get_parameters() override {
    
    model_parameters params;
    params.parameter_values = this->parameters;
    params.item_name = this->item_name;
    
    for(int i = 0; i < this->n_levels; i++){
      params.row_names.push_back(std::to_string(i));
    }
    
    for(int cl = 0; cl < this->n_classes; cl++){
      params.col_names.push_back("class_" + std::to_string(cl+1));
    }
    
    return(params);
  } 
  
};

#endif