library(rnoaa)
library(ggmap)
library(tidyverse)
library(maps)
library(mapdata)
library(mapproj)

#get list of global historical climate network daily stations
#stations<-ghcnd_stations()
write.csv(stations, "~/Dropbox/converge_diverge/datasets/LongForm/climate/ghcn_stations.csv")

#read in data
stations<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/ghcn_stations.csv")

stations_subset<-stations%>%
  select(id, latitude, longitude, name)%>%
  unique()

correloc<-read.csv("~/Dropbox/converge_diverge/datasets/LongForm/climate/siteList_LatLong.csv")

##ploting stations
map <- borders("world", colour="gray80", fill="gray80") # create a layer of borders

ggplot() +
  map+
  geom_point(data=stations_subset, aes(x=longitude, y=latitude))+
  geom_point(data=correloc, aes(x=Longitue,y=Lattitue),col="red")

#pairwise ditances
long<-nb2[,32]
lat<-nb2[,31]
longlat<-cbind(long,lat)
nb_dist1<-distm(longlat,fun=distCosine)/1000#distance btwn nb in KM
#library(geosphere)

gDistance(sp.mydata, byid=T)#library rgeos