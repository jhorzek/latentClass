#include <Rcpp.h>
#include "Distributions_Normal.hpp"

class NormalDistribution{
public:
  Normal dist;
  
  NormalDistribution(const std::vector<double> x,
                     const std::vector<std::string> par_labels,
                     const std::vector<double> par_values,
                     const std::vector<bool> free):
    dist(x,
         par_labels,
         par_values,
         free){}
  
  double get_log_likelihood(const std::vector<std::string>& par_labels,
                            const std::vector<double>& par_values){
    std::vector<double> weights(dist.data.size(), 1.0);
    return(dist.log_likelihood(par_labels,
                               par_values,
                               weights));
  }
  
  double get_log_likelihood_w(const std::vector<std::string>& par_labels,
                              const std::vector<double>& par_values,
                              const std::vector<double>& weights){
    return(dist.log_likelihood(par_labels,
                               par_values,
                               weights));
  }
  
  std::vector<double> get_gradients(const std::vector<std::string>& par_labels,
                                    const std::vector<double>& par_values){
    std::vector<double> weights(dist.data.size(), 1.0);
    return(dist.gradients(par_labels,
                          par_values,
                          weights));
  }
  
  std::vector<double> get_gradients_w(const std::vector<std::string>& par_labels,
                                      const std::vector<double>& par_values,
                                      const std::vector<double>& weights){
    return(dist.gradients(par_labels,
                          par_values,
                          weights));
  }
  
  std::vector<double> get_individual_log_likelihood(const std::vector<std::string>& par_labels,
                                                    const std::vector<double>& par_values){
    std::vector<double> weights(dist.data.size(), 1.0);
    return(dist.individual_log_likelihood(par_labels,
                                          par_values,
                                          weights));
  }
  
  std::vector<double> get_individual_log_likelihood_w(const std::vector<std::string>& par_labels,
                                                    const std::vector<double>& par_values,
                                                    const std::vector<double>& weights){
    return(dist.individual_log_likelihood(par_labels,
                                          par_values,
                                          weights));
  }
};

//RCPP_EXPOSED_CLASS(LCMR)

// Expose the LCMR class to R
RCPP_MODULE(NormalModule) {
  Rcpp::class_<NormalDistribution>("NormalDistribution")
  .constructor<std::vector<double>,
  std::vector<std::string>,
  std::vector<double>,
  std::vector<bool>>()
  .method("get_log_likelihood", &NormalDistribution::get_log_likelihood)
  .method("get_gradients", &NormalDistribution::get_gradients)
  .method("get_individual_log_likelihood", &NormalDistribution::get_individual_log_likelihood)
  
  .method("get_log_likelihood_w", &NormalDistribution::get_log_likelihood_w)
  .method("get_gradients_w", &NormalDistribution::get_gradients_w)
  .method("get_individual_log_likelihood_w", &NormalDistribution::get_individual_log_likelihood_w)
  ;
}
