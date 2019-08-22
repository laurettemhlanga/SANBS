#blood donors

rm(list = ls())

#library(PopulationSimulation)
library(inctools)
library(ggplot2)
library(glm2)





# incidence functions  ----------------------------------------------------

incidence_step <- function(donationtimes, conc = 0.02, timemin = 2001,  
                           timemax = 2011,  timedisrupt = 2007, 
                           Imin =0.01,  Ipeak =0.02, 
                           constant = T){
  
  # gives two incidence values from time before 
  # disruption and time after disruption.
  
  if (constant == T){
    
    incidence = rep(conc, length(donationtimes))
    
  }else{ 
    #varying Inicdence 
    
    incidence = ifelse(donationtimes <= timedisrupt, Imin, ifelse(donationtimes <= timemax, Ipeak, 0 ))
    
  }
  
  return(incidence)
}

linear_incidence <- function(times = 2001:2011,  timesmax = 2011, 
                             timesmin = 2001, timesdisrupt= 2009,
                             constant = 0, Imin = 0.01,
                             Idisrupt = 0.02, Ifin = 0.03)
{
  
  
  if (constant > 0) {
    
    incidence <- rep(constant,length(times))
    
  } else {
    
    incidence <-  ifelse(times < timesdisrupt, Imin + ((Idisrupt - Imin)/(timesdisrupt -timesmin)) * (times - timesmin),
                         ifelse(times <= timesmax, Idisrupt + ((Ifin - Idisrupt )/(timesmax - timesdisrupt)) * (times - timesdisrupt), 0))
    
    
  }
  
  return(incidence)
  }

# prevelence of positives functions  --------------------------------------


linear_prevalence <- function(intercept = 0.4,
                              gradient = 0,
                              times){
  
  #linear function of prevalence 
  
  prevalence <- intercept + (gradient * times)
  
  return(prevalence)
}


# recency functions  ------------------------------------------------------

linear_prevalence_recency <- function(intercept = 0.05,
                                      gradient = 0,
                                      times){
  
  #linear function of recency
  
  prevalence <- intercept + (gradient * times)
  
  return(prevalence)
}




probability_recency_disrupt <-  function(times, timesdisrupt = 2009,
                                         FRR = 0.01, MDRI = 182.5,
                                         big_T = 730,
                                         incidence = linear_incidence,
                                         prevalence = linear_prevalence
                                        
){
  
  #depends on incidence 
  
  
  MDRI = MDRI / 365; big_T = big_T / 365
  
  probabilityrecency <- FRR + (incidence(times = times, timesdisrupt = timesdisrupt ) * 
                                 ((1 - prevalence(times = times))/
                                    prevalence(times = times)) * (MDRI - (FRR * big_T)))
  
  return(probabilityrecency)
  
}



# Simulation of individual dataset  ---------------------------------------


prevalence_dataset <- function(times, 
                               prevalence_func = linear_prevalence,
                               prevalence_params = c(0.4, 0),
                               prevalences_recency_func = probability_recency_disrupt
                               # ,
                               # prevalence_recency_params = c(0.05, 0)
){
  #browser()
  #look at usage if apply functions to optimize code 
  
  #time = donationtimes
  
  prevalencedataset <- data.frame(donationtimes = times)
  
  
  prevalencedataset$HIV_status <- sapply(prevalencedataset$donationtimes, 
                                         function(x) sample(x = c(1,0), 
                                                            size = 1, prob = c(prevalence_func(times = x),
                                                                               (1 - prevalence_func(times = x))), replace = T))
  
  
  prevalencedataset$HIVrecency_status = ifelse(prevalencedataset$HIV_status == 1,
                                               sapply( prevalencedataset$donationtimes, 
                                                       function(x) sample(x = c(1,0), 
                                                                          size = 1, 
                                                                          prob = c(prevalences_recency_func(times= x),
                                                                                   (1 - prevalences_recency_func(times = x))), 
                                                                          replace = T)), NaN)
  
  
  return(prevalencedataset)
}

#prevalence_dataset()



# estimating the derivative of the incidence based on Kassanjee -----------



derivative_incidence <- function(rec_pred,
                                 prev_pred,
                                 prev_grad,
                                 prev_grad_se,
                                 rec_grad,
                                 rec_grad_se,
                                 donationtimes1,
                                 FRR, MDRI, big_T
){
  
  
  first_deriv_term <-  (((1 - 2 * (prev_pred)) / (1 - prev_pred)^2)*((rec_pred - FRR) / (MDRI - (FRR * big_T)))) * (rnorm(n = 1, mean = prev_grad, sd = prev_grad_se))
  
  
  second_deriv_term <-  (prev_pred /((1 - prev_pred) * (MDRI - (FRR * big_T)))) * (rnorm(n = 1, mean = rec_grad, sd = rec_grad_se))
  
  
  denominator <-  (1 - prev_pred) * (MDRI - FRR)
  
  return(sum(c(first_deriv_term, second_deriv_term)))
  
}



# simulations of the niterations of the donor's data set  -----------------




test_statistic_diff_inc <- function(niterations = 1,
                                    prevalence_func = linear_prevalence,
                                    prevalence_params = c(0.4, 0),
                                    prevalences_recency_func = linear_prevalence_recency,
                                    prevalence_recency_params = c(0.05, 0),
                                    start_donationtimes ,
                                    end_donationtimes,
                                    num_donations,
                                    num_interpolations,
                                    time_threshold, MDRI,
                                    big_T,  FRR)
  
{
  # browser()
  
  donationtimes <- seq(start_donationtimes,
                       end_donationtimes,
                       ((end_donationtimes - start_donationtimes)/ num_donations))

  diff_incidence_data <- data.frame(prev_1 = numeric(),
                                    prev_2 = numeric(),
                                    pre_R_1 = numeric(),
                                    pre_R_2 = numeric(),
                                    incidence_1 = numeric(),
                                    incidence_2 = numeric(),
                                    inc_rse_1 = numeric(),
                                    inc_rse_2 = numeric())

  glm_incidence_data <- data.frame(iterations = numeric(), prev_intercept = numeric(), 
                                   prev_intercept_se = numeric(),
                                   prev_gradient = numeric(), prev_gradient_se = numeric(),
                                   rec_intercept = numeric(), rec_intercept_se = numeric(),
                                   rec_gradient = numeric(), rec_gradient_se = numeric(),
                                   glm_prev = numeric(),  glm_pre_R = numeric(),
                                   glm_prev_se = numeric(), glm_pre_R_se = numeric(),
                                   glm_incidence_1 = numeric(), glm_inc_rse_1 = numeric())


  for (iterations in 1:niterations){
    
    individualdata <- prevalence_dataset(times = donationtimes,
                                         prevalence_func = prevalence_func,
                                         prevalence_params = prevalence_params,
                                         prevalences_recency_func = prevalences_recency_func
                                         # ,
                                         # prevalence_recency_params = prevalence_recency_params
    )
    
    
    # start GLM iterations  ---------------------------------------------------
    
    
    
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
    
    interpolated_times <- seq(min(donationtimes_prime), max(donationtimes_prime), (max(donationtimes_prime) - min(donationtimes_prime))/num_interpolations)
    
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
    
    
   
    
    
    glm_incidence_output <- data.frame(iterations = iterations, 
                                       prev_intercept = prev_intercept,
                                       prev_intercept_se = prev_intercept_se,
                                       prev_gradient = prev_gradient,
                                       prev_gradient_se = prev_gradient_se,
                                       rec_intercept = rec_intercept,
                                       rec_intercept_se = rec_intercept_se,
                                       rec_gradient = rec_gradient,
                                       rec_gradient_se = rec_gradient_se,
                                       glm_prev = pred_prevalences_pe, 
                                       glm_pre_R = pred_rec_prevalences_pe,
                                       glm_prev_se = pred_prevalences_se, 
                                       glm_pre_R_se = pred_rec_prevalences_se,
                                       glm_incidence = as.numeric(as.vector(glm_incidence_estimates$Incidence)), 
                                       glm_inc_rse = as.numeric(as.vector(glm_incidence_estimates$RSE)))
    
    glm_incidence_output$inci_prime_o <- derivative_incidence(rec_pred = pred_rec_prevalences_pe, prev_pre = pred_prevalences_pe,
                                                            prev_grad = prev_gradient, prev_grad_se =  prev_gradient_se, rec_grad = rec_gradient,
                                                            rec_grad_se = rec_gradient_se, FRR = FRR, MDRI = MDRI, big_T = big_T )
    
    
    
   
    
    glm_incidence_data <- rbind(glm_incidence_data, glm_incidence_output)
    
    
    # start of diff iterations --------------------------------------------------
    
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
    
    
    diff_incidence_output
    diff_incidence_data <- rbind(diff_incidence_data, diff_incidence_output )
    
    
    
    
  }
  
  diff_incidence_data$delta_inc <- diff_incidence_data$incidence_1 - diff_incidence_data$incidence_2
  diff_incidence_data$delta_inc_se <- sqrt(((diff_incidence_data$inc_rse_1  * diff_incidence_data$incidence_1)^2/niterations) +
                                             ((diff_incidence_data$inc_rse_2 * diff_incidence_data$incidence_2)^2 /niterations))
  
  diff_incidence_data$z_statistic <- diff_incidence_data$delta_inc / diff_incidence_data$delta_inc_se
  diff_incidence_data$pvalue <- 2 * pnorm(-abs(diff_incidence_data$z_statistic))
  
  glm_incidence_data$pvalue <- 2 * pnorm(-abs(glm_incidence_data$inci_prime_o / (sd(glm_incidence_data$inci_prime_o)/sqrt(niterations))))
  
  

  
  
  return(list(diff_incidence_data, glm_incidence_data))
  
}





#incidence_estimates <- rbind(incidence_estimates, incidence_estimates)

# something <- incprops()
# make data frema with (IPE_1/2, Isigma_1/2, DeltaIPE, DeltaILow, delta_I_upper, p value)
# add all these estiamtes to growing data.frame(

# 
# variosu summaries
# 
# make histogra of p values
# compre to standad normalizePath()


# ggplot(incidence_data, aes(pvalues)) +
#   geom_histogram(bins = 120)+
#   labs(x = "p-values", y = " frequency")+
#   theme_bw(base_size = 18, base_family = "")
# 
# 
# prevalence_data <- data.frame(prevalences = numeric(), prevelences_recency = numeric())
# 
# for (iterations in 1:niterations){
#   
#   individualdata <- prevalence_dataset(donationtimes = donationtimes,
#                                        prevalence_func = prevalence_func,
#                                        prevalences_recency_func = prevalence_recency
#   )
#   
#   prevalence_data <- rbind(prevalence_data, data.frame(prevalences = mean(individualdata$HIV_status), 
#                                                        prevalences_recency = (mean(individualdata$HIVrecency_status))))
#   
#   
#   
# }
# 
# point_estimate <- data.frame(mean_prev = mean( prevalence_data$prevalences),
#                              mean_prev_rece = mean(prevalence_data$prevalences_recency),
#                              sd_prev = sd( prevalence_data$prevalences),
#                              sd_prev_rece = sd(prevalence_data$prevalences_recency))
# 
# return(point_estimate)
# }
# 
# 
# 




















# 
# 
# Prevalence_data <- function(donationtimes = seq(1,365, 0.25), 
#                             prevalence_func = prevalence,
#                             prevalences_recency_func = prevalence_recency,
#                             niterations = 100
#                             ){
#   
#   prevalence_data <- data.frame(prevalences = numeric(), prevelences_recency = numeric())
#   
#   for (iterations in 1:niterations){
#     
#     individualdata <- prevalence_dataset(donationtimes = donationtimes,
#                                             prevalence_func = prevalence_func,
#                                             prevalences_recency_func = prevalence_recency
#                                             )
#     
#     prevalence_data <- rbind(prevalence_data, data.frame(prevalences = mean(individualdata$HIV_status), 
#                                                          prevalences_recency = (mean(individualdata$HIVrecency_status))))
#     
#     
#     
#   }
#   
#  point_estimate <- data.frame(mean_prev = mean( prevalence_data$prevalences),
#                               mean_prev_rece = mean(prevalence_data$prevalences_recency),
#                               sd_prev = sd( prevalence_data$prevalences),
#                               sd_prev_rece = sd(prevalence_data$prevalences_recency))
#   
#  return(point_estimate)
# }
# 
# 
# 
# dataset <- Prevalence_data()
# 
# 
# 
#   
#   
#   
# incidence <-   incprops(PrevH = dataset$mean_prev, RSE_PrevH = dataset$sd_prev, PrevR = dataset$mean_prev_rece, 
#                         RSE_PrevR = dataset$sd_prev_rece, BS_Count = 10000, Boot = F,
#                         alpha = 0.05, MDRI = 200, RSE_MDRI = 0.05, FRR = 0.001, RSE_FRR = 0,
#                         BigT = 730)
#  




