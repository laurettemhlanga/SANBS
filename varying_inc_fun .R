rm(list = ls())
source("/home/laurette/Desktop/Github/SANBS/fixing_SANBS.R", echo = FALSE)




times = seq(2001, 2011, (2011 - 2001)/10000); timesmax = 2011
timesmin = 2001; timesdisrupt= 2009; constant = 0; time_threshold = 2006
Imin = 0.002; Idisrupt = 0.003;  Ifin = 0.004
varying_con = 0.02#c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
gradient = 0; intercept = 0.5; FRR = 0.01; MDRI = 182.5
big_T = 730;  num_interpolations = 10


#c(0, 0.00025, 0.000375, 0.00075, 0.000875,0.0009375, 0.000975,0.001,0.00125, 0.00135, 0.01, 0.02, 0.03)



sensitivity_incidence <- function(times = seq(2001, 2011, (2011 - 2001)/10000), timesmax = 2011,
                        timesmin = 2001, timesdisrupt= 2009, constant = 0, time_threshold = 2006,
                        Imin = 0.02, Idisrupt = 0.03,  Ifin = 0.04,
                        varying_con = c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1),
                        gradient = 0, intercept = 0.5, FRR = 0.01, MDRI = 182.5,
                        big_T = 730,  num_interpolations = 10){

  #browser()
  
    confidenceint <-  data.frame(slope = numeric(),
                                 continoustime = numeric(),
                                 diffinc1 = numeric(),
                                 diffinc2 = numeric(),
                                 deriveinc = numeric(),
                                 derivative = numeric(), 
                                 derivative_se  =  numeric(),
                                 difference =  numeric(), 
                                 difference_se  =  numeric())
    
    if (constant > 0 ){
      
      incidencevalues <- linear_incidence(times = times, timesmax = timesmax, timesmin = timesmin,
                                                           timesdisrupt = timesdisrupt, constant = constant, Imin = Imin,
                                                           Idisrupt = Idisrupt, Ifin = Ifin, varying_con = varying_con)
    
    
      prevalence <- linear_prevalence(intercept = intercept, gradient = gradient, times = times)
      
      
      
      probabilityrecently <- probability_recency_disrupt(times = times, timesdisrupt = timesdisrupt,
                                                         FRR = FRR, MDRI = MDRI,
                                                         big_T = big_T,
                                                         incidence = incidencevalues ,
                                                         prevalence = prevalence)
      
      
      individualdata <- prevalence_dataset(times = times, 
                                           prevalence = prevalence,
                                           prevalences_recency = probabilityrecently)
      
      
      diffoutput <- diff_pvalues_calculation(individualdata = individualdata, time_threshold = time_threshold, 
                                             MDRI = MDRI, FRR = FRR, big_T = big_T)
      
      
      glmoutput <- fit_glm_derivative_calculation(individualdata = individualdata, MDRI = MDRI, FRR = FRR, big_T = big_T, 
                                                  num_interpolations = num_interpolations)
      
      
      
      
      confidenceint <-rbind(confidenceint,  data.frame(slope = 0, 
                                                       diffinc1 = diffoutput$incidence_1,
                                                       diffinc2 = diffoutput$incidence_2,
                                                       deriveinc = mean(glmoutput$glm_incidence),
                                                       derivative = mean(glmoutput$derivincidence), 
                                                       derivative_se  =  sd(glmoutput$derivincidence),
                                                       difference = diffoutput$delta_inc, 
                                                       difference_se  = diffoutput$delta_inc_se))
    
    }else{
      
    for (slope in seq_along(varying_con)){
    
    incidencevalues <- linear_incidence(times = times, timesmax = timesmax, timesmin = timesmin,
                                  timesdisrupt = timesdisrupt, constant = constant, Imin = Imin,
                                  Idisrupt = Idisrupt, Ifin = Ifin, varying_con = varying_con[slope])
    
    
    prevalence <- linear_prevalence(intercept = intercept, gradient = gradient, times = times)
    
    
    
    probabilityrecently <- probability_recency_disrupt(times = times, timesdisrupt = timesdisrupt,
                                                 FRR = FRR, MDRI = MDRI,
                                                 big_T = big_T,
                                                 incidence = incidencevalues ,
                                                 prevalence = prevalence)
    
    
    individualdata <- prevalence_dataset(times = times, 
                                   prevalence  = prevalence,
                                   prevalences_recency = probabilityrecently)
    
    
    diffoutput <- diff_pvalues_calculation(individualdata = individualdata, time_threshold = time_threshold, 
                                     MDRI = MDRI, FRR = FRR, big_T = big_T)
    
    
    glmoutput <- fit_glm_derivative_calculation(individualdata = individualdata, MDRI = MDRI, FRR = FRR, big_T = big_T, 
                                          num_interpolations = num_interpolations)
    
    
    # incidence_deriv(rec = glmoutput$glm_pre_R,
    #                 preval = glmoutput$glm_prev,
    #                 prev_gradient_se = prev_gradient_se,
    #                 prev_gradient = prev_gradient,
    #                 rec_gradient = rec_gradient,
    #                 rec_gradient_se = rec_gradient_se,
    #                 MDRI = 182.5, FRR = 0.01, big_T = 730,
    #                 prevalstderror = pred_prevalences_se,
    #                 recstderror =  pred_rec_prevalences_se)
    
    confidenceint <-rbind(confidenceint,  data.frame(slope = varying_con[slope], 
                                               diffinc1 = diffoutput$incidence_1,
                                               diffinc2 = diffoutput$incidence_2,
                                               deriveinc = mean(glmoutput$glm_incidence),
                                               derivative = mean(glmoutput$derivincidence), 
                                               derivative_se  =  sd(glmoutput$derivincidence),
                                               difference = diffoutput$delta_inc, 
                                               difference_se  = diffoutput$delta_inc_se
                                               ))
    }
    
    
    }#derivincidence 
    
    confidenceint$devlb <-  confidenceint$derivative -(1.96 * confidenceint$derivative_se)
    confidenceint$devub <-  confidenceint$derivative + (1.96 * confidenceint$derivative_se)
    
    
    confidenceint$difflb <-  confidenceint$difference -(1.96 *confidenceint$difference_se )
    confidenceint$diffub <-  confidenceint$difference + (1.96 *confidenceint$difference_se)
    
    
    return(confidenceint)

}




output3 <- data.frame()

niterations <- 1:3000

for (iterations in seq_along(niterations)){
  
iterate <- sensitivity_incidence()

output3 <- rbind(output3, iterate)
}

write.table(output3, "iterations3176")

ggplot(confidenceint,aes(x = slope,  y = derivative))+
  geom_point()+geom_errorbar(aes(ymin = devlb, ymax = devub) )


ggplot(confidenceint,aes(x = slope, y = difference))+
  geom_point()+geom_errorbar(aes(ymin = difflb, ymax = diffub))




rec <- glmoutput$glm_pre_R
preval <- glmoutput$glm_prev

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

 derivprev <-  attributes(devprev)$gradient[ ,1] *(rnorm(n = rows(gradients1), mean = prev_gradient, sd =  prev_gradient_se))
 derivrrec <-  attributes(devprev)$gradient[ ,2] * (rnorm(n = rows(gradients1), mean = rec_gradient, sd =  rec_gradient_se))
 
 derivincidence <-  derivprev +  derivrrec
 

# second derivative  ------------------------------------------------------

 
 gradpreval <- (rnorm(n = rows(gradients1), mean = prev_gradient, sd =  prev_gradient_se))
 gradrec <- (rnorm(n = rows(gradients1), mean = prev_gradient, sd =  rec_gradient_se))
 
 
 dI2dr <- deriv(~ (((rec - FRR)* gradpreval) / ((1 - preval)^2 * (MDRI - (FRR * big_T))) + 
                     (gradrec * preval)/ ((1 - preval) * (MDRI - (FRR * big_T)))), c("preval", "rec"))
 
 
 devdelta <- eval(dI2dr)
 gradients2 <- attributes(devdelta)$gradient
 
 deriv2prev <-  attributes(devdelta)$gradient[ ,1] * prevalstderror
 deriv2rrec <-  attributes(devdelta)$gradient[,2]* recstderror 
 
 incidencederiv_se <- sqrt(deriv2prev ^2 + deriv2rrec ^ 2)
 
 
 return(data.frame(derivincidence, incidencederiv_se))
  
  
}

y = incidence_deriv(rec = glmoutput$glm_pre_R,
                preval = glmoutput$glm_prev,
                prev_gradient_se = prev_gradient_se,
                prev_gradient = prev_gradient,
                rec_gradient = rec_gradient,
                rec_gradient_se = rec_gradient_se,
                MDRI = 182.5, FRR = 0.01, big_T = 730,
                prevalstderror = pred_prevalences_se,
                recstderror =  pred_rec_prevalences_se)



rec = glmoutput$glm_pre_R
preval = glmoutput$glm_prev
prev_gradient_se = prev_gradient_se
prev_gradient = prev_gradient
rec_gradient = rec_gradient
rec_gradient_se = rec_gradient_se
MDRI = 182.5; FRR = 0.01; big_T = 730
prevalstderror <- pred_prevalences_se
recstderror <-  pred_rec_prevalences_se








