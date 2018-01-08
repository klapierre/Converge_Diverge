library(tidyverse)

setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

dat<-read.csv("SpeciesRawAbundance_Oct2017.csv")%>%
  mutate(site_project_comm=paste(site_code, project_name, community_type, sep="_"))

dat1<-dat[1:90614,]
dat2<-dat[90615:196839,]
dat3<-dat[196840:267418,]

ggplot(data=dat1, aes(abundance))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")
ggplot(data=dat2, aes(abundance))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")
ggplot(data=dat3, aes(abundance))+
  geom_histogram()+
  facet_wrap(~site_project_comm, ncol=4, scales="free")

