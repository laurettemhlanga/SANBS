



diff_pvalues_calculation <- function(individualdata){
  
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
  
  return(diff_incidence_output)
}
