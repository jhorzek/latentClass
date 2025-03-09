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
  std::vector<std::string> parameter_labels;
  std::vector<double> parameter_values;
  std::vector<std::vector<std::unique_ptr<DistributionBase>>> distributions;
  
  std::vector<double> sample_weights;
  
  enum Step step = init;
  
  // For the expectation step, we need to compute the responsibilities for 
  // each class and person. These represent person-specific probabilities
  // to belong to a specific class P(C_j = c | x_i, theta_j).
  arma::mat responsibilities;
  
  void add_parameter(std::string new_parameter,
                     double start_value){
    // we only add the parameter if it is not already in the 
    // parameter vector
    for(size_t i = 0; i < this->parameter_labels.size(); i++){
      if(new_parameter.compare(parameter_labels[i]) == 0){
        return;
      }
    }
    
    parameter_labels.push_back(new_parameter);
    parameter_values.push_back(start_value);
  }
  
public:
  
  LCM(std::vector<double> class_probabilities,
      int n_persons) : class_probabilities(class_probabilities), n_persons(n_persons) {
    this->n_classes = class_probabilities.size();
    this->distributions.resize(n_classes);
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
  
  void set_parameters(std::vector<std::string> par_labels,
                      std::vector<double> new_values){
    for(size_t i = 0; i < par_labels.size(); i++){
      // fast check: if the order of the parameters is the same as before,
      // we can avoid looping over the vector every time:
      bool found = false;
      if(par_labels[i].compare(this->parameter_labels[i]) == 0){
        this->parameter_values[i] = new_values[i];
        found = true;
        continue;
      }
      // else, we need to iterate over all possible locations
      for(size_t j = 0; j < this->parameter_labels.size(); j++){
        if(par_labels[i].compare(this->parameter_labels[j]) == 0){
          this->parameter_values[j] = new_values[i];
          found = true;
          break;
        }
      }
      if(!found){
        Rcpp::stop("Could not find parameter: " + par_labels[i] + ".");
      }
    }
    this->step = parameters_changed;
  }
  
  std::vector<double> get_parameter_values(){
    return(this->parameter_values);
  }
  
  std::vector<std::string> get_parameter_labels(){
    return this->parameter_labels;
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
      for (const auto& dist : this->distributions[cl]) {
        // For each distribution, get the log-likelihood
        std::vector<double> current_ll = dist->individual_log_likelihood(this->parameter_labels, 
                                                                         this->parameter_values,
                                                                         current_weights);
        for (size_t i = 0; i < this->responsibilities.n_rows; ++i) {
          // Compute the responsibility by multiplying with the exponential of the current log-likelihood
          this->responsibilities(i, cl) *= std::exp(current_ll.at(i));
        }
      }
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
  
  double expected_log_likelihood(const std::vector<std::string>& par_labels,
                                 const std::vector<double>& par_values){
    // Compute the expected log-likelihood, which is the log-likelihood of the 
    // complete data (observed data plus latent variables).
    // Note that this likelihood is different from the likelihood of the model, which
    // would only take the observed data into account.
    if(this->step == init){
      Rcpp::stop("Please compute responsibilities once before extracting the likelihood.");
    }
    
    this->set_parameters(par_labels,
                         par_values);
    
    // We won't update the responsibilities here as those are seen as fixed
    // in the maximization step -> updating them would break this assumption.
    
    // For the log-likelihood, we must have the responsibilities (probability
    // for each person to belong to a specific class) and then have to use
    // the individual likelihoods.
    // We compute the individual log-likelihoods, while taking the new
    // class probabilities into account. Luckily we can often take a shortcut
    // so that we don't have to actually compute the likelihood for each person
    // separately.
    
    double combined_log_likelihood = 0.0;
    double resp_times_class_probs = 0.0;
    
    for(std::size_t cl = 0; cl < this->n_classes; cl++){
      
      arma::vec current_responsibilities_arma = this->responsibilities.col(cl);
      
      // Convert the Armadillo vector to a std::vector<double>
      std::vector<double> current_responsibilities(current_responsibilities_arma.begin(), 
                                                   current_responsibilities_arma.end());
      
      for(std::size_t i = 0; i < current_responsibilities.size(); i++){
        // for weighted likelihood:
        current_responsibilities.at(i) *= this->sample_weights.at(i);
      }
      
      for (const auto& dist : this->distributions[cl]) {
        for(std::size_t i = 0; i < current_responsibilities.size(); i++){
          resp_times_class_probs += current_responsibilities.at(i) *
            std::log(this->class_probabilities.at(cl));
        }
        
        combined_log_likelihood += dist->log_likelihood(this->parameter_labels,
                                                        this->parameter_values,
                                                        current_responsibilities);
      }
    }
    combined_log_likelihood += resp_times_class_probs;
    return(combined_log_likelihood);
  }
  
  // Compute the gradient of the log-likelihood
  std::vector<double> expected_log_likelihood_gradients(const std::vector<std::string>& par_labels,
                                                        const std::vector<double>& par_values){
    if(this->step == init){
      Rcpp::stop("Please compute responsibilities once before extracting the likelihood.");
    }
    
    this->set_parameters(par_labels,
                         par_values);
    
    // initialize gradients vector
    std::vector<double> grad(this->parameter_labels.size(), 0.0);
    std::vector<double> current_grad(this->parameter_labels.size(), 0.0);
    
    for(std::size_t cl = 0; cl < this->n_classes; cl++){
      
      arma::vec current_responsibilities_arma = this->responsibilities.col(cl);
      
      // Convert the Armadillo vector to a std::vector<double>
      std::vector<double> current_responsibilities(current_responsibilities_arma.begin(), 
                                                   current_responsibilities_arma.end());
      
      for(std::size_t i = 0; i < current_responsibilities.size(); i++){
        // for weighted likelihood:
        current_responsibilities.at(i) *= this->sample_weights.at(i);
      }
      
      for (const auto& dist : this->distributions[cl]) {
        
        current_grad = dist->gradients(this->parameter_labels,
                                       this->parameter_values,
                                       current_responsibilities);
        
        for(int i = 0; i < grad.size(); i++){
          grad.at(i) += current_grad.at(i);
        }
      }
    }
    
    return(grad);
    
  }
  
  // Distributions that can be added to the model
  void add_normal(int cl,
                  std::vector<double> x,
                  std::vector<std::string> dist_par_labels,
                  std::vector<double> starting_values = {0, 0}){
    this->step = init;
    if(cl > this->n_classes){
      Rcpp::stop("Class dies not exist.");
    }
    if(x.size() != this->n_persons){
      Rcpp::stop("Incorrect sample size");
    }
    std::vector<bool> free(starting_values.size());
    for(size_t i = 0; i < dist_par_labels.size(); i++){
      // skip for fixed parameters
      if(dist_par_labels[i].compare("") == 0){
        free[i] = false;
      }else{
        free[i] = true;
        add_parameter(dist_par_labels[i],
                      starting_values[i]);
      }
    }
    distributions[cl-1].push_back(std::make_unique<Normal>(x, 
                                                           dist_par_labels,
                                                           starting_values,
                                                           free));
  }
  
};

#endif