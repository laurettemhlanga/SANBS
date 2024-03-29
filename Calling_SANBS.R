
rm(list = ls())
source("/home/laurette/Desktop/Github/SANBS/SANBS.R", echo = FALSE)
#calculate I_1 - I_2

incidence_data1 <- test_statistic_diff_inc(niterations = 1,
                                            prevalence_func = linear_prevalence,
                                            prevalences_recency_func = linear_prevalence_recency,
                                            start_donationtimes = 2001,
                                            end_donationtimes = 2011,
                                            num_donations = 100000,
                                            num_interpolations = 11,
                                            time_threshold = 2006, MDRI = 182.5,
                                            big_T = 730, FRR = 0.01)


write.table(incidence_data1[[1]], file = "without_diruption_pvalues")
incidence_data1_pvalues <- read.table("/home/laurette/Desktop/Github/SANBS/without_diruption_pvalues")


write.table(incidence_data1[[2]], file = "without_diruption_glm")
incidence_data1_glm <- read.table("/home/laurette/Desktop/Github/SANBS/without_diruption_glm")

#1000 iterations 
# user      system  elapsed 
# 182.366   0.005   182.346  

#10000 iterations 

# user        system     elapsed 
# 2238.646    0.665     2238.961 

delta_inc_scatter <- ggplot(incidence_data1_glm, aes(x = iterations, y = glm_incidence)) +
                        geom_point()+
                        labs(x = "iteration", y = "Delta Incidence")+
                        # geom_hline(aes(yintercept = mean(delta_inc)), color  = "blue", linetype = "dashed", size = 1.5)+
                        # geom_hline(aes(yintercept = mean(delta_inc) - 1.96 * mean(delta_inc_se)), color  = "blue", linetype = "dashed", size = 2)+
                        # geom_hline(aes(yintercept = mean(delta_inc) + 1.96 * mean(delta_inc_se)), color  = "blue", linetype = "dashed", size = 2)+
                        theme_bw(base_size = 18, base_family = "")






inc_gradient2 <- ggplot(incidence_data1_glm, aes(x = iterations, y = inci_prime_o)) +
  geom_point()+
  labs(x = "iteration", y = "Delta Incidence")+
  geom_hline(aes(yintercept = mean(inci_prime_o)), color  = "blue", linetype = "dashed", size = 1.5)+
  geom_hline(aes(yintercept = mean(inci_prime_o) - 1.96 * (sd(inci_prime_o)/mean(inci_prime_o))), color  = "blue", linetype = "dashed", size = 2)+
  geom_hline(aes(yintercept = mean(inci_prime_o) + 1.96 * (sd(inci_prime_o)/mean(inci_prime_o))), color  = "blue", linetype = "dashed", size = 2)+
  theme_bw(base_size = 18, base_family = "")


ggplot(incidence_data1_glm, aes(pvalue)) +
  geom_histogram(bins = 10)+
  labs(x = "derivative of the incidence", y = " frequency")+
  theme_bw(base_size = 18, base_family = "")


ggplot(incidence_data1_pvalues, aes(pvalue)) +
  geom_histogram(bins = 10)+
  labs(x = "difference in incidence", y = " frequency")+
  theme_bw(base_size = 18, base_family = "")



# inc_diff2 <- ggplot(incidence_data1_pvalues, aes(x = 1:length(delta_inc), y = abs(delta_inc))) +
#   geom_point()+
#   labs(x = "iteration", y = "derivative measure")+
#   theme_bw(base_size = 18, base_family = "")

# incprops(PrevH = 0.4, RSE_PrevH = 0.0002, PrevR = 0.05, RSE_PrevR = 0.0009,
#          BS_Count = 10000, Boot = F, MDRI = 182.5, RSE_MDRI = 0.00005,
#          FRR = 0.01, RSE_FRR = 0.000002, BigT = 730)

# with disruptions  ----------------------------------------------------

incidence_data <-  test_statistic_diff_inc(niterations = 1,
                                           prevalence_func = linear_prevalence,
                                           prevalences_recency_func = probability_recency_disrupt,
                                           start_donationtimes = 2001,
                                           end_donationtimes = 2011,
                                           num_donations = 100000,
                                           num_interpolations = 10,
                                           time_threshold = 2006, MDRI = 182.5,
                                           big_T = 730, FRR = 0.01)





datan <-  data.frame(times = 2001:2011, 
                    true_incidence = linear_incidence(times  = 2001:2011, timesdisrupt = 2006),
                    timebinnning   = c(rep(incidence_data[[1]]$incidence_1, 6), rep(incidence_data[[1]]$incidence_2, 5)),
                    timebinnning_rse  = c(rep(incidence_data[[1]]$inc_rse_1, 6), rep(incidence_data[[1]]$inc_rse_2, 5)),
                    continuoustime = incidence_data[[2]]$glm_incidence,
                    continuoustime_rse = incidence_data[[2]]$glm_inc_rse)


dataq <-  data.frame(times = 2001:2011, 
                     #true_incidence = linear_incidence(times  = 2001:2011, timesdisrupt = 2006),
                     #timebinnning   = c(rep(incidence_data[[1]]$incidence_1, 6), rep(incidence_data[[1]]$incidence_2, 5)),
                     binnning  = c(rep(incidence_data[[1]]$inc_rse_1, 6), rep(incidence_data[[1]]$inc_rse_2, 5)) * 100,
                     #continuoustime = incidence_data[[2]]$glm_incidence,
                     continuous = incidence_data[[2]]$glm_inc_rse * 100)


data = melt(data = dataq, id = "times")


ggplot(data) + geom_smooth(aes(x = times, y = value, colour= variable), size  = 1, se = F) +
  scale_colour_manual(values=c("pink","purple"))+
  labs(x = "time", y = "% relative standard error", color  = "")+
  theme_bw(base_size = 18, base_family = "")



data <- data.frame(times = 2001:2011,
                   binning  = (abs(datan$true_incidence - datan$timebinnning)/ datan$true_incidence) * 100,
                   continuous = (abs(datan$true_incidence - datan$continuoustime)/ datan$true_incidence) * 100)


library(reshape2)

data = melt(data = data, id = "times")


ggplot(data) + geom_smooth(aes(x = times, y = value, colour= variable), size  = 1, se = F) +
  scale_colour_manual(values=c("pink","purple"))+
  labs(x = "time", y = "% relative bias", color  = "")+
  theme_bw(base_size = 18, base_family = "") 




write.table(incidence_data[[1]], file = "with_diruption_pvalues")
incidence_data_pvalues <- read.table("/home/laurette/Desktop/Github/SANBS/with_diruption_pvalues")


write.table(incidence_data[[2]], file = "with_diruption_glm")
incidence_data_glm <- read.table("/home/laurette/Desktop/Github/SANBS/with_diruption_glm")


# z_scores2 <- ggplot(incidence_data_pvalues, aes(z_statistic)) +
#   geom_histogram(bins = 10)+
#   labs(x = "z scores", y = " frequency")+
#   theme_bw(base_size = 18, base_family = "")



pvalues_2 <- ggplot(incidence_data_pvalues, aes(pvalue)) +
  geom_histogram(bins = 10)+
  labs(x = "p-values", y = " frequency")+
  theme_bw(base_size = 18, base_family = "")



pvalues_glm <- ggplot(incidence_data_glm, aes(pvalue)) +
  geom_histogram(bins = 10)+
  labs(x = "p-values", y = " frequency")+
  theme_bw(base_size = 18, base_family = "")


# inc_diff2 <- ggplot(incidence_data_pvalues, aes(x = 1:length(delta_inc), y = abs(delta_inc))) +
#   geom_point()+
#   labs(x = "iteration", y = "derivative measure")+
#   theme_bw(base_size = 18, base_family = "")

mean(incidence_data_glm$glm_incidence)

inc_gradient2 <- ggplot(incidence_data_glm, aes(x = iterations, y = glm_incidence)) +
  geom_point()+
  labs(x = "iteration", y = "incidence")+
  theme_bw(base_size = 18, base_family = "")




# # using linear incidence --------------------------------------------------
# 
# 
# linearincidence_data <-  test_statistic_diff_inc(niterations = 1000,
#                                            prevalence_func = linear_prevalence,
#                                            prevalences_recency_func = probability_recency_disrupt,
#                                            start_donationtimes = 2001,
#                                            end_donationtimes = 2011,
#                                            num_donations = 10000,
#                                            num_interpolations = 10,
#                                            time_threshold = 2006, MDRI = 182.5,
#                                            big_T = 730, FRR = 0.01)
# 
# 
# write.table(linearincidence_data[[1]], file = "with_diruption_pvalues_linear")
# linearincidence_data_pvalues <- read.table("/home/laurette/Desktop/Github/SANBAS/with_diruption_pvalues_linear")
# 
# 
# write.table(linearincidence_data[[2]], file = "with_diruption_glm_linear")
# linearincidence_data_glm <- read.table("/home/laurette/Desktop/Github/SANBAS/with_diruption_glm_linear")
# 
# 
# z_scores2 <- ggplot(linearincidence_data_pvalues, aes(z_statistic)) +
#   geom_histogram(bins = 10)+
#   labs(x = "z scores", y = " frequency")+
#   theme_bw(base_size = 18, base_family = "")
# 
# 
# 
# pvalues_2 <- ggplot(linearincidence_data_pvalues, aes(pvalue)) +
#   geom_histogram(bins = 10)+
#   labs(x = "p-values", y = " frequency")+
#   theme_bw(base_size = 18, base_family = "")
# 
# 
# inc_diff2 <- ggplot(linearincidence_data_pvalues, aes(x = 1:length(delta_inc), y = abs(delta_inc))) +
#   geom_point()+
#   labs(x = "iteration", y = "derivative measure")+
#   theme_bw(base_size = 18, base_family = "")
# 
# 
# 
# inc_gradient2 <- ggplot(linearincidence_data_glm[1:10,], aes(x = 1:length(inci_prime_o), y = glm_incidence)) +
#   geom_point()+
#   labs(x = "iteration", y = "derivative measure")+
#   theme_bw(base_size = 18, base_family = "")
