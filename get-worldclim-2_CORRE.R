# get worldclim data by lat long pairs
# updated July 2017 to get worldclim 2.0
# http://worldclim.org/version2

setwd("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate")

#install.packages('raster')
library(raster) 
library(data.table)
library(rgdal)
library(sp)

###tryign a different way
alberta <- as.data.frame(getData("worldclim",var="prec",res=0.5, lat=53.00124, lon=-111.519684))
plot(alberta$precip1_12)


## download bioclim 2.0 to computer ####
biotifs <- list.files("C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\wc2.0_30s_prec", full.names = T)
biolist <- lapply(biotifs, raster)
par(mfrow = c(4, 5))
lapply(biolist, image)

## for CORRE
dat <- read.csv('C:\\Users\\megha\\Dropbox\\converge_diverge\\datasets\\LongForm\\climate\\siteList_LatLong.csv')
dat
str(dat)

# intersect and extract
#sel <- dat[site_code=='nilla.au']
nutnet.pts <- SpatialPoints(dat[,.(longitude,latitude)])
layout(1)
image(biolist[[1]])
points(nutnet.pts,pch=16,cex=2)

nnextractlist <- lapply(biolist, function(x) data.table(site_code = dat$site_code, 
                                                        extract(x, nutnet.pts)))
names(nnextractlist) <- paste0('bio',1:19)
nnextract <- rbindlist(nnextractlist, idcol = 'varz')
nnbioclim2 <- dcast(nnextract, site_code ~ varz, value.var = 'V2')

nutnet.bc2 <- nnbioclim2[, c('site_code', paste0('bio',1:19)), with = F]

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

bc.names <- c('MAT','MAT_RANGE','ISO','TEMP_VAR','MAX_TEMP','MIN_TEMP','ANN_TEMP_RANGE','TEMP_WET_Q','TEMP_DRY_Q','TEMP_WARM_Q','TEMP_COLD_Q','MAP','MAP_WET_M','MAP_DRY_M','MAP_VAR','MAP_WET_Q','MAP_DRY_Q','MAP_WARM_Q','MAP_COLD_Q')

setnames(nutnet.bc2,2:20, bc.names)
nutnet.bc2

write.csv(nutnet.bc2,'~/Dropbox/NutNet-climate/worldclim2-NutNet.csv',row.names=F)

## compare to original worldclim
nn.bc <- fread('~/Dropbox/NutNet_data/climate/site-climate-table.csv')
both <- nn.bc[nutnet.bc2, on = 'site_code']
library(ggplot2)

# MAT
p <- ggplot(both, aes(x = MAT, y = i.MAT))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[MAT > 20 & i.MAT < 10] # jena.de is OFF

cat(paste0('#', bc.names), sep = '\n')

# MAT
p <- ggplot(both, aes(x = MAT, y = i.MAT))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[MAT > 20 & i.MAT < 10] # jena.de is OFF

#MAT_RANGE
p <- ggplot(both, aes(x = MAT_RANGE, y = i.MAT_RANGE))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[MAT_RANGE < 15 & i.MAT_RANGE > 15]

#ISO
p <- ggplot(both, aes(x = ISO, y = i.ISO))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)

#TEMP_VAR
p <- ggplot(both, aes(x = TEMP_VAR, y = i.TEMP_VAR/10))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[TEMP_VAR < 25 & i.TEMP_VAR > 50]

#MAX_TEMP
p <- ggplot(both, aes(x = MAX_TEMP, y = i.MAX_TEMP))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[MAX_TEMP > 30 & i.MAX_TEMP < 22]

#MIN_TEMP
#ANN_TEMP_RANGE
#TEMP_WET_Q
#TEMP_DRY_Q
#TEMP_WARM_Q
#TEMP_COLD_Q
#MAP
p <- ggplot(both, aes(x = MAP, y = i.MAP))
p + geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  geom_smooth(method = 'lm', se = F, lty = 3)
both[MAP <500 & i.MAP > 500, .(site_code, MAP, i.MAP)]

#MAP_WET_M
#MAP_DRY_M
#MAP_VAR
#MAP_WET_Q
#MAP_DRY_Q
#MAP_WARM_Q
#MAP_COLD_Q
