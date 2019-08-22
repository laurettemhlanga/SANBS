








library(parallel)




parallel_test_statistic_diff_inc <- function(niterations = 100,
                                    prevalence_func = linear_prevalence,
                                    prevalence_params = c(0.4, 0),
                                    prevalences_recency_func = linear_prevalence_recency,
                                    prevalence_recency_params = c(0.05, 0),
                                    start_donationtimes,
                                    end_donationtimes,
                                    num_donations,
                                    num_interpolations,
                                    time_threshold, MDRI ,
                                    big_T,  FRR)
  
{
  
  donationtimes <- seq(start_donationtimes,
                       end_donationtimes,
                       ((end_donationtimes - start_donationtimes)/ num_donations))
  
  numCores <- detectCores()
  
  boot_fx <- function(iterations){
    
    individualdata <- prevalence_dataset(donationtimes = donationtimes,
                                         prevalence_func = prevalence_func,
                                         prevalence_params = prevalence_params,
                                         prevalences_recency_func = prevalences_recency_func,
                                         prevalence_recency_params = prevalence_recency_params)
    
    
    glm_incidence_output <- fit_glm_derivative_calculation(individualdata)
    diff_incidence_output <- diff_pvalues_calculation(individualdata)
    
    glm_incidence_data <- rbind(data.frame(), glm_incidence_output)
    diff_incidence_data <- rbind(data.frame(), diff_incidence_output )
    
  }

   results <- mclapply(trials, boot_fx, mc.cores = numCores)
  
  diff_incidence_data$delta_inc <- diff_incidence_data$incidence_1 - diff_incidence_data$incidence_2
  diff_incidence_data$delta_inc_se <- sqrt((diff_incidence_data$inc_rse_1  * diff_incidence_data$incidence_1)^2 +
                                             (diff_incidence_data$inc_rse_2 * diff_incidence_data$incidence_2)^2)
  
  diff_incidence_data$z_statistic <- diff_incidence_data$delta_inc / diff_incidence_data$delta_inc_se
  diff_incidence_data$pvalue <- 2 * pnorm(-abs(diff_incidence_data$z_statistic))
  
  
  
  glm_incidence_data$inci_prime_o <- derivative_incidence(rec_pred = pred_rec_prevalences_pe, prev_pre = pred_prevalences_pe,
                                                          prev_grad = prev_gradient, rec_grad = rec_gradient,
                                                          FRR = FRR, MDRI = MDRI, big_T = big_T )
  
  
  return(list(diff_incidence_data, glm_incidence_data))
  
}

