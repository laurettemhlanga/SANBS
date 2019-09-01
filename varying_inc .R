rm(list = ls())
source("/home/laurette/Desktop/Github/SANBS/SANBS.R", echo = FALSE)
#calculate I_1 - I_2



linear_incidence <- function(times = 2001:2011,  timesmax = 2011, 
                             timesmin = 2001, timesdisrupt= 2009,
                             constant = 0, Imin = 0.002,
                             Idisrupt = 0.003, Ifin = 0.004, varying_con = 0.0068)
{
  
  
  if (constant > 0) {
    
    incidence <- rep(constant,length(times))
    
  } else {
    
    incidence <-  ifelse(times < timesdisrupt, Imin + (((Idisrupt + varying_con) - Imin)/(timesdisrupt -timesmin)) * (times - timesmin),
                         ifelse(times <= timesmax, (Idisrupt + varying_con) + ((Ifin - (Idisrupt  + varying_con))/(timesmax - timesdisrupt)) * (times - timesdisrupt), 0))
    
    
  }
  
  return(incidence)
}








incidence_data <-   test_statistic_diff_inc(niterations = 10,
                                            prevalence_func = linear_prevalence,
                                            prevalences_recency_func = probability_recency_disrupt,
                                            start_donationtimes = 2001,
                                            end_donationtimes = 2011,
                                            num_donations = 100000,
                                            num_interpolations = 10,
                                            time_threshold = 2006, MDRI = 182.5,
                                            big_T = 730, FRR = 0.01)


mean(incidence_data[[1]]$delta_inc)
mean(incidence_data[[1]]$delta_inc_se)

mean(incidence_data[[2]]$inci_prime_o)
sd(incidence_data[[2]]$inci_prime_o)



confidenceint <-  data.frame(slope = c(0, 0.00025, 0.000375, 0.00075, 0.000875,0.0009375, 0.000975,0.001,0.00125, 0.00135),
                             derivative = c(3.623505e-05, 1.010127e-05,1.228036e-05, 1.26231e-05, 1.609608e-05, 8.343776e-06 ,1.473744e-05,1.728447e-05,5.777134e-05, 4.477765e-05), 
                             derivative_se  = c( 2.016921e-05, 8.272812e-06, 1.021584e-05, 1.369102e-05, 1.022329e-05, 1.058574e-05,1.12579e-05 ,1.8521e-05,2.141365e-05, 1.535327e-05   ),
                             difference = c(-0.00218,-0.001663, -0.001203, -0.002435, -0.003155, -0.002461, -0.00321,-0.00367,-0.00759 ,-0.007341), 
                             difference_se  = c( 0.001833668,0.001485135, 0.001552504, 0.001569387, 0.001641596, 0.001652844, 0.001657524,0.001669743,0.002014730,   0.002064964))


confidenceint$devlb <-  confidenceint$derivative -(1.96 *confidenceint$derivative_se)
confidenceint$devub <-  confidenceint$derivative + (1.96 *confidenceint$derivative_se)


confidenceint$difflb <-  confidenceint$difference -(1.96 *confidenceint$difference_se )
confidenceint$diffub <-  confidenceint$difference + (1.96 *confidenceint$difference_se)


ggplot(confidenceint,aes(x = slope, 
                         y = derivative))+
  geom_point()+geom_errorbar(aes(ymin = devlb, ymax = devub) )


ggplot(confidenceint,aes(x = slope, 
                         y = difference))+
  geom_point()+geom_errorbar(aes(ymin = difflb, ymax = diffub) )



derivative_delta_method <- function(predprevalence, predrecency,
                                    MDRI, FRR , bigT, prevstderror,
                                    recstderror, devprevalence,
                                    devrecency
){
  
  # MDRI = MDRI / 365; bigT = bigT / 365
  
  partialdevprevalence <- ((((predrecency - FRR) * devprevalence  *(2 * predprevalence))/ ((MDRI - (FRR * bigT))*(1 - predprevalence)^3)) - 
    ((devrecency * predprevalence) / ((MDRI - (FRR * bigT))*(1 - predprevalence)^2)))^2
  
  
  partialdevrecency <- (((1 - (2 * predprevalence)) * devprevalence) / ((MDRI - (FRR * bigT))*(1 - predprevalence)^2))^2
  
  
  overall_error <-  sqrt(partialdevprevalence * prevstderror^2 + partialdevrecency * recstderror ^2)
  
  return(overall_error)
  
}



derivative_delta_method(predprevalence = 0.4016502, predrecency = 0.05039398,
                        MDRI = 182.5, FRR = 0.01, bigT = 730, prevstderror = 0.003098620,
                        recstderror = 0.002188974, devprevalence = -0.0005708268,
                        devrecency = 6.680063e-05)



























# # probability of recency  -------------------------------------------------
# 
# varying_con = c(0.001,0.002, 0.003, 0.005, 0.006)
# 
# Incidence <- sapply(varying_con, function(x) linear_incidence(times = 2001:2011,  timesmax = 2011, 
#                                     timesmin = 2001, timesdisrupt= 2009,
#                                     constant = 0, Imin = 0.01,
#                                     Idisrupt = 0.02, Ifin = 0.03, varying_con = x))
# 
# 
# 
# 
# 
# 
# 
# 
# probability_recency_disrupt(times = 2001:2011, timesdisrupt = 2009,
#                             FRR = 0.01, MDRI = 182.5,
#                             big_T = 730, incidence = linear_incidence,
#                             prevalence = linear_prevalence)



# probability_recency_disrupt(times = 2001:2011, timesdisrupt = 2009,
#                             FRR = 0.01, MDRI = 182.5,
#                             big_T = 730,
#                             incidence = varing_incidence,
#                             prevalence = linear_prevalence)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# incgradient <- c(0, 0.001, 0.005, 0.01, 0.02)
# 
# probability_recency <- list()
# 
# #varyslope <- 1
# for (varyslope in  seq_along(incgradient) ){
# 
# 
# # incidence <- varing_incidence(times = 2001:2011, min_times = 2001, max_times = 2011,
# #                                max_inc = 0.05, min_inc = 0.01,
# #                                varying_con = incgradient[varyslope],
# #                                intercept = 0.01)
# 
# 
# probability_recency <- probability_recency_disrupt(times = times, timesdisrupt = 2009,
#                                                    FRR = 0.01, MDRI = 182.5,
#                                                    big_T = 730,
#                                                    incidence = varing_incidence(times, min_times = 2001, max_times = 2011,
#                                                                                 max_inc = 0.05, min_inc = 0.01,
#                                                                                 varying_con = incgradient[varyslope],
#                                                                                 intercept = 0.01),
#                                                    prevalence = linear_prevalence)
#   
# }
