#ifndef MODEL_H
#define MODEL_H
#include "Distributions_Normal.hpp"
#include "RcppArmadillo.h"

inline std::vector<double> elementwise_product(const std::vector<double>& a,
                                               const std::vector<double>& b){
  if(a.size() != b.size())
    Rcpp::stop("Vectors must be of same size.");
  std::vector prod(a.size(), 0.0);
  for(int i = 0; i < a.size(); i++){
    prod.at(i) = a.at(i)*b.at(i);
  }
  return(prod);
}

enum Step {
  init,
  parameters_changed,
  responsibilities_updated,
  class_probabilities_updated
};

class LCM {
  
private:
  int n_classes;
  int n_persons;
  std::vector<double> class_probabilities;
  std::vector<std::unique_ptr<DistributionBase>> distributions;
  
  std::vector<double> sample_weights;
  
  enum Step step = init;
  
  // For the expectation step, we need to compute the responsibilities for 
  // each class and person. These represent person-specific probabilities
  // to belong to a specific class P(C_j = c | x_i, theta_j).
  arma::mat responsibilities;
  
public:
  
  LCM(std::vector<double> class_probabilities,
      int n_persons) : class_probabilities(class_probabilities), n_persons(n_persons) {
    this->n_classes = class_probabilities.size();
    this->responsibilities.resize(n_persons, n_classes);
    this->step = init;
    this->sample_weights = std::vector<double>(n_persons, 1.0);
  }
  
  void set_sample_weights(const std::vector<double> sample_weights){
    if(sample_weights.size() != this->n_persons){
      Rcpp::stop("sample_weights must have the length n_persons.");
    }
    this->sample_weights = sample_weights;
  }
  
  int get_n_classes(){
    return(this->n_classes);
  }
  
  int get_n_persons(){
    return(this->n_persons);
  }
  
  void set_class_probabilities(std::vector<double> new_class_probabilities){
    if(new_class_probabilities.size() != this->n_classes){
      std::string cl = std::to_string(this->n_classes);
      Rcpp::stop("new_class_probabilities must be of length " + cl + ".");
    }
    this->class_probabilities = new_class_probabilities;
    this->step = parameters_changed;
  }
  
  std::vector<double> get_class_probabilities(){
    return(this->class_probabilities);
  }
  
  
  model_parameters get_parameters(){
    model_parameters params;
    for(const auto& dist: this->distributions){
      model_parameters current_pars = dist->get_parameters();
      params.parameter_names.insert(params.parameter_names.end(), 
                                    current_pars.parameter_names.begin(), 
                                    current_pars.parameter_names.end());
      params.parameter_values.insert(params.parameter_values.end(), 
                                     current_pars.parameter_values.begin(), 
                                     current_pars.parameter_values.end());
    }
    return(params);
  }
  
  void update_responsibilities(){
    
    this->responsibilities.fill(0.0);
    
    // responsibilities are independent of the weights
    std::vector<double> current_weights(n_persons, 1.0);
    
    for (int cl = 0; cl < this->n_classes; cl++) {
      // Add class probabilities to responsibilities
      for (size_t i = 0; i < this->responsibilities.n_rows; ++i) {
        this->responsibilities(i, cl) = this->class_probabilities.at(cl);
      }
    }
    // Now, we need to iterate over all distributions and add the likelihoods
    // returned by the distributions.
    for (const auto& dist : this->distributions) {
      // For each distribution, get the log-likelihood matrix.
      // This matrix has the individuals in the rows and the classes in the
      // columns.
      arma::mat current_ll = dist->individual_log_likelihood(current_weights);
      // We now have to multiply those with our responsibilities element-wise to
      // update the responsibilities:
      this->responsibilities = this->responsibilities % arma::exp(current_ll);
    }
    
    // the responsibilities can become extremely small, which can result in division
    // by zero. To prevent this from happening, we will specify a lower bound:
    for (size_t i = 0; i < this->responsibilities.n_rows; ++i) {
      for (size_t j = 0; i < this->responsibilities.n_cols; ++i) {
        if (this->responsibilities(i,j) < 1e-8) { // Avoid division by zero
          this->responsibilities(i,j) = 1e-8;
        }
      }
    }
    
    // Normalize all of the responsibilities
    arma::vec row_sums = arma::sum(this->responsibilities, 1);
    for (std::size_t i = 0; i < this->responsibilities.n_rows; ++i) {
      if (row_sums(i) > 0) { // Avoid division by zero
        this->responsibilities.row(i) /= row_sums(i);
      } else {
        Rcpp::stop("Division by zero.");
      }
    }
    
    this->step = responsibilities_updated;
  }
  
  void update_class_probabilities(){
    // Sum the responsibilities for every class
    // class probabilities have to use the weights
    std::vector<double> class_probs(this->responsibilities.n_cols,
                                    0.0);
    for(std::size_t cl = 0; cl < this->responsibilities.n_cols; cl++){
      double n = 0.0;
      for(std::size_t i = 0; i < this->n_persons; i++){
        n += this->sample_weights.at(i);
        class_probs.at(cl) += this->sample_weights.at(i)*this->responsibilities(i, cl);
      }
      class_probs.at(cl) /= n;
    }
    this->class_probabilities = class_probs;
  }
  
  arma::mat get_responsibilities(){
    if(this->step == init){
      Rcpp::stop("Please compute responsibilities once before extracting the responsibilities.");
    }
    return(this->responsibilities);
  }
  
  /*
   double log_likelihood(const std::vector<std::string>& par_labels,
   const std::vector<double>& par_values){
   if(this->step == init){
   Rcpp::stop("Please compute responsibilities once before extracting the likelihood.");
   }
   
   this->set_parameters(par_labels,
   par_values);
   
   std::vector<double> ind_lik(this->n_persons, 0.0);
   
   for(std::size_t cl = 0; cl < this->n_classes; cl++){
   
   std::vector<double> ind_cl_lik(this->n_persons, this->class_probabilities.at(cl));
   
   for (const auto& dist : this->distributions[cl]) {
   std::vector<double> ind_cl_dist_lik = dist->individual_log_likelihood(this->parameter_labels,
   this->parameter_values,
   this->sample_weights);
   for(int i = 0; i < this->n_persons; i++){
   // we currently have the log-likelihood, so we first have to compute the 
   // likelihood:
   ind_cl_dist_lik.at(i) = std::exp(ind_cl_dist_lik.at(i));
   }
   ind_cl_lik = elementwise_product(ind_cl_lik,
   ind_cl_dist_lik);
   }
   
   for(int i = 0; i < this->n_persons; i++){
   ind_lik.at(i) += ind_cl_lik.at(i);
   }
   }
   
   // Now, we can finally sum everything up
   double ll = 0.0;
   for(int i = 0; i < this->n_persons; i++){
   ll += std::log(ind_lik.at(i));
   }
   
   return(ll);
   }
   */
  
  void maximize(){
    // Executes the maximization step for each of the distributions based 
    // on the complete data (observed data plus latent variables).
    if(this->step == init){
      Rcpp::stop("Please compute responsibilities once before extracting the likelihood.");
    }
    
    // We won't update the responsibilities here as those are seen as fixed
    // in the maximization step -> updating them would break this assumption.
    
    for(const auto& dist: this->distributions){
      dist->maximize_parameters(this->responsibilities);
    }
  }
  
  // Distributions that can be added to the model
  void add_normal(std::string item_name,
                  std::vector<double> x,
                  std::vector<double> initial_means,
                  std::vector<double> initial_sds,
                  bool sd_equal){
    if(x.size() != this->n_persons)
      Rcpp::stop("Expected size of x to the be same as the number of persons.");
    if(initial_means.size() != this->n_classes)
      Rcpp::stop("Expected one mean per class.");
    if(initial_sds.size() != this->n_classes)
      Rcpp::stop("Expected one sd per class.");
    
    this->step = init;
    
    this->distributions.push_back(std::make_unique<Normal>(item_name,
                                                           x, 
                                                           initial_means,
                                                           initial_sds,
                                                           sd_equal));
  }
  
};

#endif