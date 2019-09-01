#blood donors

rm(list = ls())

#library(PopulationSimulation)
library(inctools)
library(ggplot2)
library(glm2)





# incidence functions  ----------------------------------------------------


linear_incidence <- function(times,  timesmax = 2011, 
                             timesmin = 2001, timesdisrupt,
                             constant = 0, Imin = 0.002,
                             Idisrupt = 0.003, Ifin = 0.004, 
                             varying_con = 0.02)
{
  
  
  if (constant > 0) {
    
    incidence <- rep(constant,length(times))
    
  } else {
    
    incidence <-  ifelse(times < timesdisrupt, Imin + (((Idisrupt + varying_con) - Imin)/(timesdisrupt -timesmin)) * (times - timesmin), 
                         Imin + (((Idisrupt + varying_con) - Imin)/(timesdisrupt -timesmin)) * (timesdisrupt -timesmin))
    
    
  }
  
  return(incidence)
}

#incidence = linear_incidence()

  
  # prevelence of positives functions  --------------------------------------


linear_prevalence <- function(intercept,
                              gradient,
                              times ){
  
  #linear function of prevalence 
  
  prevalence <- intercept + (gradient * times)
  
  return(prevalence)
}

#prevalence = linear_prevalence()


probability_recency_disrupt <-  function(times, timesdisrupt,
                                         FRR, MDRI,
                                         big_T,
                                         incidence,
                                         prevalence 
                                         
){
  
  #depends on incidence 
  
  
  MDRI = MDRI / 365; big_T = big_T / 365
  
  probabilityrecency <- FRR + (incidence * ((1 - prevalence) / prevalence) * (MDRI - (FRR * big_T)))
  
  return(probabilityrecency)
  
}


#probrecency <- probability_recency_disrupt()


# Simulation of individual dataset  ---------------------------------------


prevalence_dataset <- function(times, 
                               prevalence,
                               prevalences_recency 
                             
){
  
  prevalencedataset <- data.frame(donationtimes = times)
  
  
  prevalencedataset$HIV_status <- sapply(seq_along(prevalence), function(x) sample(x = c(1,0), 
                                                                                   size = 1, prob = c(prevalence[x], (1 - prevalence[x])),  replace = T))
  
  
  prevalencedataset$HIVrecency_status = ifelse(prevalencedataset$HIV_status == 1,
                                               sapply(seq_along(prevalences_recency), 
                                                      function(x) sample(x = c(1,0), size = 1, 
                                                                         prob = c(prevalences_recency[x], (1 - prevalences_recency[x])), replace = T)), NaN)
  
  return(prevalencedataset)
}

#individual <- prevalence_dataset()

derivative_incidence <- function(rec_pred,
                                 prev_pred,
                                 prev_grad,
                                 prev_grad_se,
                                 rec_grad,
                                 rec_grad_se,
                                 donationtimes1,
                                 FRR, MDRI, big_T
){
  
  #MDRI = MDRI / 365; big_T = big_T / 365
  
  first_deriv_term <- ((rec_pred - FRR) * (rnorm(n = length(rec_pred), mean = prev_grad, sd = prev_grad_se)))/((1 - prev_pred)^2 *(MDRI - (FRR * big_T)))
  
  
  second_deriv_term <- (prev_pred * (rnorm(n = length(rec_pred), mean = rec_grad, sd = rec_grad_se))) / ((1 - prev_pred) * (MDRI - (FRR * big_T)))
  
  derivative = first_deriv_term + second_deriv_term
  
  return(derivative)
  
}


# error term derivative  --------------------------------------------------


derivative_delta_method <- function(predprevalence, predrecency,
                                    MDRI, FRR , bigT, prevstderror,
                                    recstderror, devprevalence,
                                    devrecency
){
  
  #MDRI = MDRI / 365; bigT = bigT / 365
  
  partialdevprevalence <- ((((predrecency - FRR) * devprevalence  * 2)/ ((MDRI - (FRR * bigT))*(1 - predprevalence)^3)) +
                             (devrecency  / ((MDRI - (FRR * bigT))*(1 - predprevalence)^2)))
  
  
  partialdevrecency <- (devprevalence / ((MDRI - (FRR * bigT)) * (1 - predprevalence)^2)) 
  
  
  overall_error <-  sqrt((partialdevprevalence * prevstderror)^2 + (partialdevrecency * recstderror)^2)
  
  return(overall_error)
  
}


# inbuilt functions -------------------------------------------------------

incidence_deriv <- function(rec, preval,
                            prev_gradient_se,
                            prev_gradient,
                            rec_gradient,
                            rec_gradient_se,
                            MDRI, FRR, big_T,
                            prevalstderror,
                            recstderror
){
  
  dIdpdr <-  deriv(~ ((preval * (rec - FRR)) / ((1 - preval)* (MDRI - (FRR *big_T)))), c("preval", "rec"))
  
  
  devprev <- eval(dIdpdr)
  
  gradients1 <- attributes(devprev)$gradient
  
  derivprev <-  attributes(devprev)$gradient[ ,1] *(rnorm(n = nrow(gradients1), mean = prev_gradient, sd =  prev_gradient_se))
  derivrrec <-  attributes(devprev)$gradient[ ,2] * (rnorm(n = nrow(gradients1), mean = rec_gradient, sd =  rec_gradient_se))
  
  derivincidence <-  derivprev +  derivrrec
  
  
  # second derivative  ------------------------------------------------------
  
  
  gradpreval <- (rnorm(n = nrow(gradients1), mean = prev_gradient, sd =  prev_gradient_se))
  gradrec <- (rnorm(n = nrow(gradients1), mean = prev_gradient, sd =  rec_gradient_se))
  
  
  dI2dr <- deriv(~ (((rec - FRR)* gradpreval) / ((1 - preval)^2 * (MDRI - (FRR * big_T))) + 
                      (gradrec * preval)/ ((1 - preval) * (MDRI - (FRR * big_T)))), c("preval", "rec"))
  
  
  devdelta <- eval(dI2dr)
  gradients2 <- attributes(devdelta)$gradient
  
  deriv2prev <-  attributes(devdelta)$gradient[ ,1] * prevalstderror
  deriv2rrec <-  attributes(devdelta)$gradient[,2]* recstderror 
  
  incidencederiv_se <- sqrt(deriv2prev ^2 + deriv2rrec ^ 2)
  
  
  return(data.frame(derivincidence, incidencederiv_se))
  
  
}

# y = incidence_deriv(rec = glmoutput$glm_pre_R,
#                     preval = glmoutput$glm_prev,
#                     prev_gradient_se = prev_gradient_se,
#                     prev_gradient = prev_gradient,
#                     rec_gradient = rec_gradient,
#                     rec_gradient_se = rec_gradient_se,
#                     MDRI = 182.5, FRR = 0.01, big_T = 730,
#                     prevalstderror = pred_prevalences_se,
#                     recstderror =  pred_rec_prevalences_se)




# estimating the derivative of the incidence based on Kassanjee -----------


diff_pvalues_calculation <- function(individualdata, 
                                     time_threshold, 
                                     MDRI, FRR,
                                     big_T){
  
  data1 <- individualdata[individualdata$donationtimes <= time_threshold, ] # subsett
  data2 <- individualdata[individualdata$donationtimes > time_threshold,] # subsett
  
  prev_1 <- mean(data1$HIV_status)
  sigma_prev_1 <- sqrt(( prev_1 *(1 -  prev_1 ))/nrow(data1))
  pre_R_1 <- sum(data1$HIVrecency_status, na.rm = T)/sum(data1$HIV_status) #dealing with missing values 
  sigma_prev_R_1 <-  sqrt((  pre_R_1 *(1 - pre_R_1 ))/ nrow(data1[which(data1$HIV_status == 1),]))
  
  prev_2 <- mean(data2$HIV_status)
  sigma_prev_2 <- sqrt(( prev_2 *(1 -  prev_2 ))/nrow(data2))
  pre_R_2 <-  sum(data2$HIVrecency_status, na.rm = T)/sum(data2$HIV_status)
  sigma_prev_R_2 <- sqrt((  pre_R_2 *(1 - pre_R_2 ))/ nrow(data2[which(data1$HIV_status == 1),]))
  
  
  diff_incidence_estimates = incprops(PrevH = c(prev_1, prev_2), RSE_PrevH = c(sigma_prev_1 / prev_1, sigma_prev_2 / prev_2),
                                      PrevR = c(pre_R_1, pre_R_2), RSE_PrevR = c(sigma_prev_R_1 / pre_R_1 , sigma_prev_R_2 / pre_R_2), 
                                      BS_Count = 10000, Boot = F, alpha = 0.05, MDRI = MDRI, RSE_MDRI = 0, 
                                      FRR = FRR, RSE_FRR = 0, BigT = big_T)$Incidence.Statistics
  
  
  diff_incidence_output <- data.frame(prev_1 = prev_1,
                                      prev_2 = prev_2,
                                      prev_R_1 = pre_R_1,
                                      prev_R_2 = pre_R_2,
                                      incidence_1 = as.numeric(as.vector(diff_incidence_estimates$Incidence[1])), 
                                      incidence_2 = as.numeric(as.vector(diff_incidence_estimates$Incidence[2])), 
                                      inc_rse_1 = as.numeric(as.vector(diff_incidence_estimates$RSE[1])),
                                      inc_rse_2 = as.numeric(as.vector(diff_incidence_estimates$RSE[2])))
  
  
  
  diff_incidence_output$delta_inc <- diff_incidence_output$incidence_1 - diff_incidence_output$incidence_2
  diff_incidence_output$delta_inc_se <- sqrt(((diff_incidence_output$inc_rse_1  * diff_incidence_output$incidence_1)^2) +
                                             ((diff_incidence_output$inc_rse_2 * diff_incidence_output$incidence_2)^2))
  
  diff_incidence_output$z_statistic <- diff_incidence_output$delta_inc / diff_incidence_output$delta_inc_se
  diff_incidence_output$pvalue <- 2 * pnorm(-abs(diff_incidence_output$z_statistic))
  
  return(diff_incidence_output)
}

#diff_pvalues_calculation()

# estimating the derivative of the incidence based on Kassanjee -----------


fit_glm_derivative_calculation <- function(individualdata,  
                                           MDRI, FRR, 
                                           big_T, 
                                           num_interpolations
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
                                     glm_incidence = as.numeric(as.vector(glm_incidence_estimates$Incidence)), 
                                     glm_inc_rse = as.numeric(as.vector(glm_incidence_estimates$RSE)))
  
  
  
  
  # glm_incidence_output$inci_prime_o <- derivative_incidence(rec_pred = pred_rec_prevalences_pe, prev_pre = pred_prevalences_pe,
  #                                                           prev_grad = prev_gradient, prev_grad_se =  prev_gradient_se, rec_grad = rec_gradient,
  #                                                           rec_grad_se = rec_gradient_se, FRR = FRR, MDRI = MDRI, big_T = big_T)
  # 
  # 
  # 
  # 
  # glm_incidence_output$inci_prime_o_se <- derivative_delta_method(predprevalence= pred_prevalences_pe, predrecency = pred_rec_prevalences_pe,
  #                                                                 MDRI = MDRI, FRR = FRR, bigT = big_T, prevstderror =  pred_prevalences_se,
  #                                                                 recstderror = pred_rec_prevalences_se, devprevalence = prev_gradient,
  #                                                                 devrecency = rec_gradient)
  
  
  output <- incidence_deriv(rec =   glm_incidence_output$glm_pre_R,
                  preval =   glm_incidence_output$glm_prev,
                  prev_gradient_se = prev_gradient_se,
                  prev_gradient = prev_gradient,
                  rec_gradient = rec_gradient,
                  rec_gradient_se = rec_gradient_se,
                  MDRI = MDRI, FRR = FRR, big_T = big_T,
                  prevalstderror = pred_prevalences_se,
                  recstderror =  pred_rec_prevalences_se)
  
  
  glm_incidence_output <- cbind(glm_incidence_output, output)

  
  glm_incidence_output$pvalue <- 2 * pnorm(-abs( glm_incidence_output$derivincidence/ sd(glm_incidence_output$derivincidence)))
  
  
  
  
  return(glm_incidence_output)
}



#fit_glm_derivative_calculation()










# derivative_delta_method(predprevalence = 0.4016502, predrecency = 0.05039398,
#                         MDRI = 182.5, FRR = 0.01, bigT = 730, prevstderror = 0.003098620,
#                         recstderror = 0.002188974, devprevalence = -0.0005708268,
#                         devrecency = 6.680063e-05)

# simulations of the niterations of the donor's data set  -----------------














