


fit_glm_derivative_calculation <- function(individualdata
                                          ){
  
  individualdata$donationtimes_prime <- individualdata$donationtimes - ((max(individualdata$donationtimes)+min(individualdata$donationtimes))/2)
  
  prevelance_fit <- glm2(formula = HIV_status ~ 1 + donationtimes_prime,
                         family = binomial(link = "identity"), data = individualdata)
  
  recency_fit <- glm2(formula = HIVrecency_status ~ 1 + donationtimes_prime,
                      family = binomial(link = "identity"), 
                      data = individualdata[individualdata$HIV_status == 1, ])
  
  prev_intercept <- summary(prevelance_fit)$coefficients[1,1] 
  prev_gradient <- summary(prevelance_fit)$coefficients[2,1] 
  
  prev_intercept_se <- summary(prevelance_fit)$coefficients[1,2] 
  prev_gradient_se <- summary(prevelance_fit)$coefficients[2,2]
  
  
  donationtimes_prime <- individualdata$donationtimes_prime
  
  interpolated_times <- seq(min(donationtimes_prime), 
                            max(donationtimes_prime), 
                            (max(donationtimes_prime - min(donationtimes_prime)))/num_interpolations)
  
  prev_predictions <- predict(prevelance_fit, newdata = data.frame(donationtimes_prime = interpolated_times),
                              se.fit = TRUE)
  pred_prevalences_pe <- prev_predictions$fit
  pred_prevalences_se <- prev_predictions$se.fit
  
  
  
  
  rec_intercept <- summary(recency_fit)$coefficients[1,1] 
  rec_gradient <- summary(recency_fit)$coefficients[2,1] 
  
  rec_intercept_se <- summary(recency_fit)$coefficients[1,2] 
  rec_gradient_se <- summary(recency_fit)$coefficients[2,2]
  
  
  rec_predictions <- predict(recency_fit, newdata = data.frame(donationtimes_prime = interpolated_times), 
                             se.fit = TRUE)
  pred_rec_prevalences_pe <- rec_predictions$fit
  pred_rec_prevalences_se <- rec_predictions$se.fit
  
  
  
  glm_incidence_estimates = incprops(PrevH = pred_prevalences_pe, RSE_PrevH = (pred_prevalences_se / pred_prevalences_pe),
                                     PrevR = pred_rec_prevalences_pe, RSE_PrevR = (pred_rec_prevalences_se / pred_rec_prevalences_pe), 
                                     BS_Count = 10000, Boot = F, alpha = 0.05, MDRI = MDRI, RSE_MDRI = 0, 
                                     FRR = FRR, RSE_FRR = 0, BigT = big_T)$Incidence.Statistics
  
  
  glm_incidence_output <- data.frame(glm_prev = pred_prevalences_pe, 
                                     glm_pre_R = pred_rec_prevalences_pe,
                                     glm_prev_se = pred_prevalences_se, 
                                     glm_pre_R_se = pred_rec_prevalences_se,
                                     glm_incidence = as.numeric(as.vector(glm_incidence_estimates$Incidence[1])), 
                                     glm_inc_rse = as.numeric(as.vector(glm_incidence_estimates$RSE[1])))
  
  
  return(glm_incidence_output)
}
