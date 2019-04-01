# get worldclim data by lat long pairs
# updated July 2017 to get worldclim 2.0
# http://worldclim.org/version2

setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm")

#install.packages('raster')
library(raster) 
library(data.table)
library(rgdal)
library(sp)
library(tidyverse)

###tryign a different way
r <- getData("worldclim", var="bio",res=2.5)

r <- r[[c(1,12)]] #bio1 and bio12 are MAT and MAP
names(r) <- c("Temp","Prec")

dat <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\siteList_LatLong.csv')


coords <- dat[,3:4]
points <- SpatialPoints(coords, proj4string=r@crs)

values <- extract(r,pts)

df <- cbind.data.frame(coordinates(points),values)



## field definitions from http://www.worldclim.org/bioclim

# BIO1 = Annual Mean Temperature
# BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
# BIO3 = Isothermality (BIO2/BIO7) (* 100)
# BIO4 = Temperature Seasonality (standard deviation *100)
# BIO5 = Max Temperature of Warmest Month
# BIO6 = Min Temperature of Coldest Month
# BIO7 = Temperature Annual Range (BIO5-BIO6)
# BIO8 = Mean Temperature of Wettest Quarter
# BIO9 = Mean Temperature of Driest Quarter
# BIO10 = Mean Temperature of Warmest Quarter
# BIO11 = Mean Temperature of Coldest Quarter
# BIO12 = Annual Precipitation
# BIO13 = Precipitation of Wettest Month
# BIO14 = Precipitation of Driest Month
# BIO15 = Precipitation Seasonality (Coefficient of Variation)
# BIO16 = Precipitation of Wettest Quarter
# BIO17 = Precipitation of Driest Quarter
# BIO18 = Precipitation of Warmest Quarter
# BIO19 = Precipitation of Coldest Quarter