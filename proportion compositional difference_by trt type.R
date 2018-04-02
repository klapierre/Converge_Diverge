library(tidyverse)

setwd('C:\\Users\\lapie\\Dropbox (Smithsonian)\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

prop <- read.csv('proportions significant comp difference_by trt type.csv')%>%
  #calculate proportions from data
  mutate(proportion=comp_change/total_possible)%>%
  #calculate expected proportions
  mutate(expected_change=(total_possible*sum(comp_change))/sum(total_possible), expected_change_prop=expected_change/total_possible,
         expected_nochange=(total_possible*sum(no_change))/sum(total_possible), expected_nochange_prop=expected_nochange/total_possible)%>%
  #calculate chi-squared
  mutate(chi=(comp_change-expected_change)^2/expected_change)

#calculate chi-squared
sum(prop$chi) #14.27226, df=19
pchisq(14.27226, df=19)

chisq.test(x = prop$comp_change,
           p = prop$expected_change_prop)

