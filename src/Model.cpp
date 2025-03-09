#include <RcppArmadillo.h>
#include "Model.hpp"
// [[Rcpp::depends(RcppArmadillo)]]

// Latent Class Model
class LCMR: LCM{
public:
  
  LCMR(std::vector<double> class_probabilities,
       int n_persons):
  LCM(class_probabilities, n_persons){}
  
  std::vector<double> get_class_probabilities_R(){
    return(this->get_class_probabilities());
  }
  
  int get_n_classes_R(){
    return(this->get_n_classes());
  }
  
  int get_n_persons_R(){
    return(this->get_n_persons());
  }
  
  Rcpp::NumericVector get_parameters_R(){
    
    std::vector<double> pars_values = this->get_parameter_values();
    
    std::vector<std::string> pars_names = this->get_parameter_labels();
    
    if (pars_values.size() != pars_names.size()) {
      Rcpp::stop("pars_values and pars_names must have the same size.");
    }
    
    Rcpp::NumericVector pars_vector(pars_values.begin(), pars_values.end());
    
    // Set the names of the vector
    pars_vector.attr("names") = Rcpp::StringVector(pars_names.begin(), pars_names.end());
    
    return pars_vector;
  }
  
  void set_class_probabilities_R(std::vector<double> new_class_probabilities){
    return(this->set_class_probabilities(new_class_probabilities));
  }
  
  void set_parameters_R(std::vector<std::string> parameter_names,
                        std::vector<double> parameter_values){
    this->set_parameters(parameter_names,
                         parameter_values);
  }
  
  void update_responsibilities_R(){
    this->update_responsibilities();
  }
  
  void update_class_probabilities_R(){
    this->update_class_probabilities();
  }
  
  arma::mat get_responsibilities_R(){
    return(this->get_responsibilities());
  }
  
  double log_likelihood_R(std::vector<std::string> parameter_names,
                          std::vector<double> parameter_values){
    return(this->log_likelihood(parameter_names,
                                parameter_values));
  }
  
  Rcpp::NumericVector gradients_R(std::vector<std::string> parameter_names,
                                  std::vector<double> parameter_values){
    
    std::vector<double> grad = this->gradients(parameter_names,
                                               parameter_values);
    
    std::vector<std::string> grad_names = this->get_parameter_labels();
    
    if (grad.size() != grad_names.size()) {
      Rcpp::stop("grad and grad_names must have the same size.");
    }
    
    Rcpp::NumericVector grad_vector(grad.begin(), grad.end());
    
    // Set the names of the vector
    grad_vector.attr("names") = Rcpp::StringVector(grad_names.begin(), grad_names.end());
    
    return grad_vector;
  }
  
  void add_normal_R(int cl,
                    std::vector<double> data,
                    Rcpp::StringVector parameter_names,
                    std::vector<double> parameter_values = {0.0, 0.0}){
    if(parameter_names.length() != 2){
      Rcpp::stop("parameter_names must be of length 2");
    }
    if(parameter_values.size() != 2){
      Rcpp::stop("starting_values must be of size 2.");
    }
    
    // Names are currently in a StringVector, but we need them as a std::vector<std::string>
    // for the model.
    std::vector<std::string> parameter_names_vec(parameter_names.size());
    for (int i = 0; i < parameter_names.size(); i++){
      if(Rcpp::StringVector::is_na(parameter_names(i))){
        parameter_names_vec[i] = "";
      }else{
        parameter_names_vec[i] = Rcpp::as< std::string >(parameter_names(i));
      }
    }
    this->add_normal(cl, data, parameter_names_vec, parameter_values);
  }
  
  
  
};

//RCPP_EXPOSED_CLASS(LCMR)

// Expose the LCMR class to R
RCPP_MODULE(LCMModule) {
  Rcpp::class_<LCMR>("LCMR")
  .constructor<std::vector<double>,int>()
  .method("update_class_probabilities", &LCMR::update_class_probabilities_R)
  .method("get_class_probabilities", &LCMR::get_class_probabilities_R)
  .method("set_class_probabilities", &LCMR::set_class_probabilities_R)
  .method("get_n_classes", &LCMR::get_n_classes_R)
  .method("get_n_persons", &LCMR::get_n_persons_R)
  .method("get_parameters", &LCMR::get_parameters_R)
  .method("set_parameters", &LCMR::set_parameters_R)
  .method("update_responsibilities", &LCMR::update_responsibilities_R)
  .method("get_responsibilities", &LCMR::get_responsibilities_R)
  .method("log_likelihood", &LCMR::log_likelihood_R)
  .method("gradients", &LCMR::gradients_R)
  .method("add_normal", &LCMR::add_normal_R);
}
