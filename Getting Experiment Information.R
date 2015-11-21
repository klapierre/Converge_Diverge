setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\datasets\\LongForm')

#setwd("~/Dropbox/converge_diverge/datasets/LongForm")

library(tidyr)
library(dplyr)
library(vegan)

siteInfo <- read.csv("SiteInfo_11202015.csv")%>%
  select(-X, -species_num)
species <- read.csv("SpeciesRawAbundance_11202015.csv")
