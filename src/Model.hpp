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
  init,  // status while model options are changed
  ready, // status after the responsibilities have been initialized and we are ready to go
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
  
  // functions for expectation maximization
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
    
    this->step = ready;
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
  
  void expectation(){
    this->update_responsibilities();
  }
  
  void maximization(){
    // Executes the maximization step for each of the distributions based 
    // on the complete data (observed data plus latent variables).
    if(this->step == init){
      this->update_responsibilities();
    }
    
    // We won't update the responsibilities here as those are seen as fixed
    // in the maximization step -> updating them would break this assumption.
    
    for(const auto& dist: this->distributions){
      dist->maximize_parameters(this->responsibilities);
    }
  }
  
  
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
    this->step = init;
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
    this->step = init;
  }
  
  std::vector<double> get_class_probabilities(){
    return(this->class_probabilities);
  }
  
  std::vector<model_parameters> get_parameters(){
    std::vector<model_parameters> params;
    for(const auto& dist: this->distributions){
      params.push_back(dist->get_parameters());
    }
    return(params);
  }
  
  arma::mat get_responsibilities(){
    if(this->step == init){
      this->update_responsibilities();
    }
    return(this->responsibilities);
  }
  
  double log_likelihood(){
    if(this->step == init){
      this->update_responsibilities();
    }
    
    std::vector<double> ind_lik(this->n_persons, 0.0);
    arma::mat ind_likelihood(this->n_persons, this->n_classes);
    ind_likelihood.fill(1.0);
    
    std::vector<double> identity_weights(this->n_persons, 1.0);
    for (const auto& dist : this->distributions) {
      arma::mat ind_class_lik = arma::exp(dist->individual_log_likelihood(identity_weights));
      // element wise product to get the individual likelihoods:
      ind_likelihood = ind_likelihood % ind_class_lik;
    }
    
    // take class probabilities and sample weights into account
    for(int i = 0; i < this->n_persons; i++){
      for(int cl = 0; cl < this->n_classes; cl++){
        ind_likelihood(i,cl) = ind_likelihood(i,cl) * this->sample_weights.at(i) * this->class_probabilities.at(cl);
      }
    }
    
    arma::colvec ind_log_likelihood = arma::sum(ind_likelihood, 1);
    // Now, we can finally sum everything up
    double ll = 0.0;
    for(int i = 0; i < this->n_persons; i++){
      ll += std::log(ind_log_likelihood(i));
    }
    
    return(ll);
  }
  
  void expectation_maximization(int max_iter = 1000,
                                double convergence_criterion =1e-7){
    // initialize responsibilities
    this->update_responsibilities();
    double ll_old = this->log_likelihood();
    double ll_new = 0.0;
    int n_iter = 0;
    bool converged = false;
    
    for(n_iter = 0; n_iter < max_iter; n_iter++){
      // Expectation
      this->expectation();
      if(n_iter > 0)
        this->update_class_probabilities();
      // Maximization
      this->maximization();
      ll_new = this->log_likelihood();
      
      if(std::abs((ll_new - ll_old)/ll_new) < convergence_criterion){
        break;
      }
      ll_old = ll_new;
    }
    
    if(n_iter == max_iter){
      Rcpp::warning("EM algorithm did not converge.");
    }else{
      converged = true;
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
    
    this->distributions.push_back(std::make_unique<Normal>(item_name,
                                                           x, 
                                                           initial_means,
                                                           initial_sds,
                                                           sd_equal));    
    this->step = init;
  }
  
};

#endif