
Varinginc <- read.table("/home/laurette/Desktop/Github/SANBS/iterations3176")


head(Varinginc)

varying_slope <-  c(0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)

y = lapply(varying_slope, function(x) filter(Varinginc, slope == x))


data <- data.frame(slope = numeric(), 
                   diffmean = numeric(),
                   diffstderror = numeric(),
                   dermean = numeric(), 
                   derstderror = numeric())



for(len in 1:11){
  
  
  
  data <- rbind(data, data.frame( slope =  mean(y[[len]]$slope),
                                  diffmean = mean(y[[len]]$difference),
                                  diffstderror = mean(y[[len]]$difference_se), 
                                  dermean = mean(y[[len]]$derivative), 
                                  derstderror = mean(y[[len]]$derivative_se)))
  
}


summarydata <-  data 
summarydata$difflb <- summarydata$diffmean - (1.96 * summarydata$diffstderror)
summarydata$diffub <- summarydata$diffmean + (1.96 * summarydata$diffstderror)

summarydata$derlb <- summarydata$dermean - (1.96 * summarydata$derstderror)
summarydata$derub <- summarydata$dermean + (1.96 * summarydata$derstderror)

summarydata

write.csv(summarydata, " results3176")















summarydata <-  data 
summarydata$difflb <- summarydata$diffmean - (1.96 * summarydata$diffstderror)
summarydata$diffub <- summarydata$diffmean + (1.96 * summarydata$diffstderror)

summarydata$derlb <- summarydata$dermean - (1.96 * summarydata$derstderror)
summarydata$derub <- summarydata$dermean + (1.96 * summarydata$derstderror)

summarydata

write.csv(summarydata, " results3176")
summarytable <- read.csv("/Users/laurette/Desktop/Github/SANBS/ results3176")

ggplot(summarytable,aes(x = slope,  y =  dermean))+
  geom_point()+geom_errorbar(aes(ymin = derlb, ymax = derub))+
  labs(x = "Incidence Slope", y = " Estimated incidence slope ")+
  theme_bw(base_size = 18, base_family = "")



ggplot(summarytable,aes(x = slope, y = diffmean))+
  geom_point()+geom_errorbar(aes(ymin = difflb, ymax = diffub))






ggplot(summarytable,aes(x = slope,  y =  dermean))+
  geom_point()+geom_errorbar(aes(ymin = derlb, ymax = derub))+
  labs(x = "Incidence Slope", y = " Estimated incidence slope ")+
  theme_bw(base_size = 18, base_family = "")



















