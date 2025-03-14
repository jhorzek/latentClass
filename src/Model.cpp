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
  
  Rcpp::List get_parameters_R(){
    
    std::vector<model_parameters> params = this->get_parameters();
    
    Rcpp::List par_list = Rcpp::List::create();
    
    for(int i = 0; i < params.size(); i++){
      Rcpp::NumericMatrix current_values(params.at(i).parameter_values.n_rows, 
                                         params.at(i).parameter_values.n_cols);
      Rcpp::CharacterVector current_row_names;
      Rcpp::CharacterVector current_col_names;
      for(int r = 0; r < params.at(i).parameter_values.n_rows; r++){
        current_row_names.push_back(params.at(i).row_names.at(r));
        for(int c = 0; c < params.at(i).parameter_values.n_cols; c++){
          if(r == 0){
            current_col_names.push_back(params.at(i).col_names.at(c));
          }
          current_values(r,c) = params.at(i).parameter_values(r,c);
        }
      }
      
      Rcpp::rownames(current_values) = current_row_names;
      Rcpp::colnames(current_values) = current_col_names;
      
      par_list.push_back(current_values,
                         params.at(i).item_name);
    }
    
    return par_list;
  }
  
  void set_class_probabilities_R(std::vector<double> new_class_probabilities){
    return(this->set_class_probabilities(new_class_probabilities));
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
  
  void maximize_R(){
    this->maximize();
  }
  
  double log_likelihood_R(){
    return(this->log_likelihood());
  }
  
  void add_normal_R(std::string item_name,
                    std::vector<double> x,
                    std::vector<double> initial_means,
                    std::vector<double> initial_sds,
                    bool sd_equal){
    
    this->add_normal(item_name, x, initial_means, initial_sds, sd_equal);
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
  .method("update_responsibilities", &LCMR::update_responsibilities_R)
  .method("get_responsibilities", &LCMR::get_responsibilities_R)
  .method("add_normal", &LCMR::add_normal_R)
  .method("maximize", &LCMR::maximize_R)
  .method("log_likelihood", &LCMR::log_likelihood_R);
}
