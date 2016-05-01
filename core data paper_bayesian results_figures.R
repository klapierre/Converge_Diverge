library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)

setwd('C:\\Users\\Kim\\Dropbox\\working groups\\converge diverge working group\\converge_diverge\\converge diverge core paper_2015\\bayesian output')

theme_set(theme_bw())
theme_update(axis.title.x=element_text(size=20, vjust=-0.35, margin=margin(t=15)), axis.text.x=element_text(size=16),
             axis.title.y=element_text(size=20, angle=90, vjust=0.5, margin=margin(r=15)), axis.text.y=element_text(size=16),
             plot.title = element_text(size=24, vjust=2),
             panel.grid.major=element_blank(), panel.grid.minor=element_blank(),
             legend.title=element_blank(), legend.text=element_text(size=20))


###bar graph summary statistics function
#barGraphStats(data=, variable="", byFactorNames=c(""))

barGraphStats <- function(data, variable, byFactorNames) {
  count <- length(byFactorNames)
  N <- aggregate(data[[variable]], data[byFactorNames], FUN=length)
  names(N)[1:count] <- byFactorNames
  names(N) <- sub("^x$", "N", names(N))
  mean <- aggregate(data[[variable]], data[byFactorNames], FUN=mean)
  names(mean)[1:count] <- byFactorNames
  names(mean) <- sub("^x$", "mean", names(mean))
  sd <- aggregate(data[[variable]], data[byFactorNames], FUN=sd)
  names(sd)[1:count] <- byFactorNames
  names(sd) <- sub("^x$", "sd", names(sd))
  preSummaryStats <- merge(N, mean, by=byFactorNames)
  finalSummaryStats <- merge(preSummaryStats, sd, by=byFactorNames)
  finalSummaryStats$se <- finalSummaryStats$sd / sqrt(finalSummaryStats$N)
  return(finalSummaryStats)
}  

##################################################################################
##################################################################################

mean <- read.csv('mean change_experiment_coefs.csv')%>%
  mutate(yr10=Intercepts + 10*Slopes + (10^2)*Quads)

ggplot(data=barGraphStats(data=mean, variable='yr10', byFactorNames=c('plot_mani')), aes(x=plot_mani, y=mean, fill=plot_mani)) +
  geom_bar(stat='identity') +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se, width=0.2)) +
  scale_fill_manual(values=c('')) #finish this

mean <- read.csv('mean change_experiment_coefs.csv')%>%
  mutate(curve1='stat_function(fun=function(x){',
         curve2=' + ',
         curve3='*x + ',
         curve4='*x^2}, xlim=c(0,',
         curve5='), colour=',
         curve6=') +',
         color=ifelse(plot_mani==1, '#1400E599', ifelse(plot_mani==2, '#4A06AC99', ifelse(plot_mani==3, '#800C7499', ifelse(plot_mani==4, '#B6123C99', '#EC180499')))),
         curve=paste(curve1, Intercepts, curve2, Slopes, curve3, Quads, curve4, experiment_length, curve5, color, curve6))

# mean2 <- mean%>%
#   select(Quads, Slopes, Intercepts, color, experiment_length)

# curves <- function(data, intercept, slope, quad, color, experiment_length) {
#   for (i in 1:nrow(data)) {
#   p <- p + stat_function(fun=function(x){intercept + slope*x + quad*x^2}, colour=color, xlim=c(0,experiment_length))
#   }
#   return(p)
# }

# p <- ggplot(data=data.frame(x=c(0,0))) +
#   ylim(c(0,1)) + xlim(c(0,28))
# curveOutput <- curves(mean2, Intercepts, Slopes, Quads, color, experiment_length)
  
# meanLines <- alply(as.matrix(mean), 1, paste(mean$curve1, mean$Intercepts, mean$curve2, mean$Slopes, mean$curve3, mean$Quads, mean$curve4, mean$experiment_length, mean$curve5, mean$color, mean$curve6))

# meanLines <- alply(as.matrix(mean2), 1, function(dat) { 
#   stat_function(fun=function(x) {dat['Quads'] * x^2 + dat['Slopes'] * x + dat['Intercepts']})  
# })
# 
# ggplot(data=data.frame(x=c(0,0))) +
#   meanLines +
#   scale_y_continuous(limits=c(0,1)) +
#   scale_x_continuous(limits=c(0,28))

p <- ggplot(data=data.frame(x=c(0,0))) +
  ylim(c(-1,2)) + xlim(c(0,28))

for (i in 1:nrow(mean)) {
  p <- p + stat_function(aes(y=0), fun=function(x) {mean[i,'Quads']*x^2 + mean[i,'Slopes']*x + mean[i,'Intercepts']}, xlim=c(0, mean[i,'experiment_length']), colour='grey')
}


for (i in 1:nrow(mean)) {
  print(mean[i,'curve'], quote=FALSE,max.levels=0)
}

print(p)



p <- ggplot(data=data.frame(x=c(0,0))) +
  coord_cartesian(xlim=c(0,28), ylim=c(0,1)) +
  xlim(0,28) + ylim(-1,2)

p <- p + 
  stat_function(fun=function(x){ 0.230929386  +  0.021502518 *x +  -0.001364054 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.262171061  +  0.02401176 *x +  -0.001525823 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.264315375  +  0.030791937 *x +  -0.001286028 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.206609096  +  0.031686224 *x +  -0.001175413 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.227275001  +  0.021500391 *x +  -0.00082508 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.270091532  +  0.027955621 *x +  -0.001133344 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.248493062  +  0.029207605 *x +  -0.001166054 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.25755068  +  0.029005406 *x +  -0.001172893 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.285460706  +  0.034425311 *x +  -0.001183539 *x^2}, size=1, xlim=c(0, 13 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.252928172  +  0.001981279 *x +  -0.000406796 *x^2}, size=1, xlim=c(0, 13 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.260585623  +  0.04300529 *x +  -0.001497189 *x^2}, size=1, xlim=c(0, 14 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.219919105  +  0.051229869 *x +  -0.0015955 *x^2}, size=1, xlim=c(0, 11 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.37634329  +  0.056351027 *x +  -0.002208445 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.433885054  +  0.065284875 *x +  -0.002671376 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.530320075  +  0.054212157 *x +  -0.003370371 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.505870741  +  0.009981243 *x +  -0.001773763 *x^2}, size=1, xlim=c(0, 10 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.488987074  +  0.013030401 *x +  -0.001834057 *x^2}, size=1, xlim=c(0, 10 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.484917572  +  0.019967268 *x +  -0.001608783 *x^2}, size=1, xlim=c(0, 10 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.387016459  +  0.039793413 *x +  -0.002103177 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.41234162  +  0.035825734 *x +  -0.001983006 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.389915554  +  0.035089413 *x +  -0.001925995 *x^2}, size=1, xlim=c(0, 8 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.370785749  +  0.024994854 *x +  -0.001472578 *x^2}, size=1, xlim=c(0, 8 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.310387372  +  0.042219887 *x +  -0.001975358 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.413064501  +  0.045393826 *x +  -0.00241283 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.375913542  +  0.036219329 *x +  -0.001965967 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.338939008  +  0.032617514 *x +  -0.001255627 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.287855876  +  0.02973926 *x +  -0.001099468 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.332686233  +  0.031268726 *x +  -0.001218994 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.350133114  +  0.030973806 *x +  -0.00123076 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.351048582  +  0.037485093 *x +  -0.001483422 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.176098128  +  0.023011823 *x +  -0.001517118 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.197645645  +  0.016056216 *x +  -0.001328941 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.139780634  +  0.018570418 *x +  -0.00128736 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.178500188  +  0.034384225 *x +  -0.001931856 *x^2}, size=1, xlim=c(0, 16 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.21031503  +  0.04433613 *x +  -0.002585511 *x^2}, size=1, xlim=c(0, 16 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.168717393  +  0.023587201 *x +  -0.001663977 *x^2}, size=1, xlim=c(0, 16 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.171530146  +  0.012670687 *x +  -0.001056369 *x^2}, size=1, xlim=c(0, 16 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.056500464  +  0.03138637 *x +  -0.001115602 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.036331593  +  0.021532631 *x +  -0.000853358 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.043494051  +  0.016490855 *x +  -0.000756066 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.036586787  +  0.017710923 *x +  -0.000888481 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.058007813  +  0.025308658 *x +  -0.001014379 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.031204735  +  0.028832469 *x +  -0.000926348 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.205895487  +  0.011412366 *x +  -0.002044338 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.228201187  +  0.008703459 *x +  -0.001977204 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.324029763  +  0.029693835 *x +  -0.00294958 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.297351986  +  0.02500285 *x +  -0.002734868 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.285805899  +  0.019490479 *x +  -0.002490953 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.21825915  +  0.009963729 *x +  -0.002228367 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.233517739  +  0.029103656 *x +  -0.002951409 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.209729359  +  0.015962118 *x +  -0.002418104 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.229964751  +  0.015391477 *x +  -0.002460041 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.250137875  +  0.0079407 *x +  -0.00226389 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.283103219  +  0.01690313 *x +  -0.002331721 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.256035165  +  0.00629114 *x +  -0.001919981 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.280592275  +  0.025900795 *x +  -0.002667701 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.300746635  +  0.031266174 *x +  -0.002934133 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.240956573  +  0.014305837 *x +  -0.002156138 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.246521514  +  0.031779115 *x +  -0.00242456 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.209583558  +  0.009936593 *x +  -0.001562067 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.273619681  +  0.029979143 *x +  -0.002424654 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.298308124  +  0.041813192 *x +  -0.00289906 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.200072034  +  0.017587664 *x +  -0.001831191 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.199834363  +  0.01524691 *x +  -0.001158153 *x^2}, size=1, xlim=c(0, 14 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.193679582  +  0.021933711 *x +  -0.00151169 *x^2}, size=1, xlim=c(0, 14 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.092756796  +  0.00947595 *x +  -0.000472847 *x^2}, size=1, xlim=c(0, 14 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.300178112  +  0.004406443 *x +  -0.000272609 *x^2}, size=1, xlim=c(0, 23 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.414275994  +  0.054686563 *x +  -0.001778748 *x^2}, size=1, xlim=c(0, 23 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.162265736  +  0.007515903 *x +  -0.000300267 *x^2}, size=1, xlim=c(0, 23 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.238162534  +  0.090281176 *x +  -0.002742005 *x^2}, size=1, xlim=c(0, 23 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.23110513  +  -0.001098398 *x +  0.000290017 *x^2}, size=1, xlim=c(0, 23 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.431619097  +  0.034674722 *x +  -0.001197686 *x^2}, size=1, xlim=c(0, 23 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.332504534  +  0.022716138 *x +  -0.001433731 *x^2}, size=1, xlim=c(0, 23 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.526499813  +  0.051668293 *x +  -0.002729824 *x^2}, size=1, xlim=c(0, 23 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.212155578  +  0.009646862 *x +  -0.000636513 *x^2}, size=1, xlim=c(0, 10 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.297137079  +  0.042198876 *x +  -0.001606892 *x^2}, size=1, xlim=c(0, 10 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.274940653  +  0.00633864 *x +  -0.001127782 *x^2}, size=1, xlim=c(0, 10 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.488647665  +  0.092501368 *x +  -0.004392884 *x^2}, size=1, xlim=c(0, 10 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.224272382  +  0.010444898 *x +  -0.000534399 *x^2}, size=1, xlim=c(0, 10 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.354865101  +  0.079669083 *x +  -0.00290311 *x^2}, size=1, xlim=c(0, 10 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.185333106  +  0.018662883 *x +  -0.001116384 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.216815344  +  0.023667588 *x +  -0.001387654 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.18501105  +  0.043916505 *x +  -0.001759745 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.134870692  +  0.029399766 *x +  -0.001019035 *x^2}, size=1, xlim=c(0, 4 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.084538683  +  0.010773588 *x +  -0.001416068 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.103025822  +  0.02063256 *x +  -0.001783486 *x^2}, size=1, xlim=c(0, 6 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.289375134  +  0.033612517 *x +  -0.000861912 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.243143675  +  0.030087941 *x +  -0.000573828 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.239141204  +  0.05581129 *x +  -0.001529445 *x^2}, size=1, xlim=c(0, 8 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.170299003  +  0.038571371 *x +  -0.001190205 *x^2}, size=1, xlim=c(0, 8 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.231644356  +  0.045839225 *x +  -0.00168189 *x^2}, size=1, xlim=c(0, 8 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.300760687  +  0.037862866 *x +  -0.001476611 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.330861924  +  0.033613708 *x +  -0.001375621 *x^2}, size=1, xlim=c(0, 6 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.183128566  +  0.006690919 *x +  -0.000390798 *x^2}, size=1, xlim=c(0, 8 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.221359678  +  0.031098681 *x +  -0.001162498 *x^2}, size=1, xlim=c(0, 8 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.206548019  +  0.015010719 *x +  -0.00032087 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.177079861  +  0.012995269 *x +  -5.44e-05 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.168546157  +  0.010486112 *x +  2.88e-05 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.160470783  +  0.022087083 *x +  -0.000720115 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.111038079  +  0.019262931 *x +  -0.000642474 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.141435145  +  0.02610461 *x +  -0.000973803 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.126222151  +  0.020588836 *x +  -0.000727862 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.122275143  +  0.019326985 *x +  -0.000683863 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.148797591  +  0.028556604 *x +  -0.001083467 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.213510916  +  0.032284485 *x +  -0.001890864 *x^2}, size=1, xlim=c(0, 13 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.245896798  +  0.023463602 *x +  -0.001633979 *x^2}, size=1, xlim=c(0, 13 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.239492862  +  0.029539103 *x +  -0.001843644 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.24326889  +  0.021690511 *x +  -0.001472961 *x^2}, size=1, xlim=c(0, 13 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.249007619  +  0.025506519 *x +  -0.00170458 *x^2}, size=1, xlim=c(0, 13 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.236549194  +  0.035860355 *x +  -0.002274106 *x^2}, size=1, xlim=c(0, 13 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.225421778  +  0.032975041 *x +  -0.00194032 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.251787793  +  0.033789962 *x +  -0.002161701 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.225248501  +  0.033361235 *x +  -0.001924413 *x^2}, size=1, xlim=c(0, 13 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.2261655  +  0.027389684 *x +  -0.001734446 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.243098165  +  0.030633377 *x +  -0.001811821 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.203956434  +  0.022703548 *x +  -0.001422709 *x^2}, size=1, xlim=c(0, 13 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.213013031  +  0.019043464 *x +  -0.001372079 *x^2}, size=1, xlim=c(0, 13 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.206904101  +  0.02252959 *x +  -0.001641347 *x^2}, size=1, xlim=c(0, 13 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.387045892  +  0.03220648 *x +  -0.001141764 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.351790432  +  0.032924428 *x +  -0.001095016 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.351972216  +  0.026820346 *x +  -0.000793766 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.337911425  +  0.02847979 *x +  -0.000839261 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.357430231  +  0.02758678 *x +  -0.00088651 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.360967481  +  0.02104657 *x +  -0.000604602 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.377678382  +  0.027230103 *x +  -0.000919029 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.344820816  +  0.034229322 *x +  -0.001125056 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.298235368  +  -0.000747876 *x +  0.000243963 *x^2}, size=1, xlim=c(0, 24 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.646820613  +  0.024231738 *x +  -0.000792657 *x^2}, size=1, xlim=c(0, 24 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.171654001  +  0.005861324 *x +  -0.001018593 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.244530742  +  0.018780528 *x +  -0.001734863 *x^2}, size=1, xlim=c(0, 7 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.196132001  +  0.023560746 *x +  -0.000659286 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.160551169  +  0.013483499 *x +  -0.000241208 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.402715514  +  0.018596363 *x +  -0.000976014 *x^2}, size=1, xlim=c(0, 25 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.345973614  +  0.021489928 *x +  -0.000880253 *x^2}, size=1, xlim=c(0, 25 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.356586047  +  0.033082393 *x +  -0.001035054 *x^2}, size=1, xlim=c(0, 25 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.373853782  +  0.03420627 *x +  -0.001352192 *x^2}, size=1, xlim=c(0, 25 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.378325215  +  0.041800687 *x +  -0.001948674 *x^2}, size=1, xlim=c(0, 25 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.346667233  +  0.021443568 *x +  -0.000586184 *x^2}, size=1, xlim=c(0, 25 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.386339088  +  0.013074915 *x +  -0.000804871 *x^2}, size=1, xlim=c(0, 25 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.452625602  +  0.022434998 *x +  -0.001418926 *x^2}, size=1, xlim=c(0, 25 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.445655219  +  0.033188043 *x +  -0.002082931 *x^2}, size=1, xlim=c(0, 25 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.361731279  +  0.021366924 *x +  -0.001116561 *x^2}, size=1, xlim=c(0, 25 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.337818449  +  0.010623636 *x +  -0.000257857 *x^2}, size=1, xlim=c(0, 25 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.301573124  +  0.010723715 *x +  -0.000383004 *x^2}, size=1, xlim=c(0, 25 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.193902622  +  0.029048534 *x +  -0.001235522 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.18373831  +  0.02216798 *x +  -0.000930215 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.20742189  +  0.020333642 *x +  -0.001039395 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.176285781  +  0.013252862 *x +  -0.000687049 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.149509414  +  0.02430702 *x +  -0.001013057 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.131342203  +  0.015690021 *x +  -0.000666863 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.134903017  +  0.025878508 *x +  -0.000777135 *x^2}, size=1, xlim=c(0, 19 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.167058713  +  0.004461523 *x +  0.000105696 *x^2}, size=1, xlim=c(0, 19 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.186414875  +  0.006639574 *x +  -0.000399876 *x^2}, size=1, xlim=c(0, 12 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.216915636  +  0.036364619 *x +  -0.001814552 *x^2}, size=1, xlim=c(0, 12 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.249747172  +  0.044406393 *x +  -0.002199385 *x^2}, size=1, xlim=c(0, 12 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.175268831  +  0.005226502 *x +  -0.000395803 *x^2}, size=1, xlim=c(0, 16 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.142193073  +  0.012028609 *x +  0.000137172 *x^2}, size=1, xlim=c(0, 16 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.161395778  +  0.028078369 *x +  -0.000343627 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.143497116  +  0.03288198 *x +  -0.000527793 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.20518477  +  0.00015646 *x +  0.00155233 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.241345234  +  0.032013213 *x +  0.000496587 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.31587413  +  0.031379055 *x +  0.000320753 *x^2}, size=1, xlim=c(0, 11 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.314182991  +  0.003520255 *x +  0.000291995 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.32631611  +  0.031664277 *x +  -0.000511443 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.359924331  +  0.045986217 *x +  -0.000865623 *x^2}, size=1, xlim=c(0, 11 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.266884163  +  0.056284421 *x +  -0.001675281 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.280226497  +  0.053399277 *x +  -0.001549594 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.279816399  +  0.036711939 *x +  -0.001209437 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.296853353  +  0.038498471 *x +  -0.00123178 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.23399189  +  0.032650008 *x +  -0.001264682 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.167682834  +  0.02130891 *x +  -0.000699754 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.177982399  +  0.025512814 *x +  -0.000869724 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.181873519  +  0.025245795 *x +  -0.000866869 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.246126728  +  0.029025311 *x +  -0.001080994 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.250579107  +  0.016959069 *x +  -0.000719076 *x^2}, size=1, xlim=c(0, 3 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.165482207  +  0.016323005 *x +  -0.000525092 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.251976637  +  0.01974712 *x +  -0.000780308 *x^2}, size=1, xlim=c(0, 9 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.215662461  +  0.016847421 *x +  -0.000734618 *x^2}, size=1, xlim=c(0, 9 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.228364766  +  0.014081113 *x +  -0.00065569 *x^2}, size=1, xlim=c(0, 9 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.243304873  +  0.043823614 *x +  -0.000778101 *x^2}, size=1, xlim=c(0, 9 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.253866777  +  0.042833119 *x +  -0.000772412 *x^2}, size=1, xlim=c(0, 9 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.217004273  +  0.039462425 *x +  -0.000530521 *x^2}, size=1, xlim=c(0, 9 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.179172881  +  0.007500438 *x +  0.000490577 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.176647902  +  0.033122629 *x +  -0.000405396 *x^2}, size=1, xlim=c(0, 9 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.17276163  +  0.009962308 *x +  6.57e-05 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.203789707  +  0.021451819 *x +  -0.000184531 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.1528175  +  0.003900707 *x +  -0.000262252 *x^2}, size=1, xlim=c(0, 9 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.197601496  +  0.017488478 *x +  -0.000720524 *x^2}, size=1, xlim=c(0, 9 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.283507447  +  0.026167983 *x +  -0.001276188 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.290673479  +  0.016324596 *x +  -0.000258838 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.334444182  +  0.02048846 *x +  -0.000486821 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.338435254  +  0.023714714 *x +  -0.00062812 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.256341823  +  0.031811269 *x +  -0.000669143 *x^2}, size=1, xlim=c(0, 5 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.299758658  +  0.029523195 *x +  -0.00065442 *x^2}, size=1, xlim=c(0, 5 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.252979381  +  0.038472101 *x +  -0.000924434 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.154122394  +  0.027309213 *x +  -0.000447056 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.209291785  +  0.043498156 *x +  -0.001206248 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.206741457  +  0.028264916 *x +  -0.000591282 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.208998907  +  0.025900114 *x +  -0.001134547 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.336826168  +  0.026762927 *x +  -0.001416737 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.226807314  +  0.018823826 *x +  -0.000883024 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.294033029  +  0.037523202 *x +  -0.001579189 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.237639555  +  0.014729341 *x +  -0.000608712 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.259423382  +  0.0258683 *x +  -0.001655734 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.359927819  +  0.008296517 *x +  -0.001248443 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.432920844  +  0.044692722 *x +  -0.002807794 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.381273052  +  0.034273045 *x +  -0.001820549 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.315438784  +  0.019653038 *x +  -0.001562072 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.345452932  +  0.02497614 *x +  -0.001950477 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.356919927  +  0.016319271 *x +  -0.001616396 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.303416411  +  0.018223354 *x +  -0.001606864 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.235814574  +  0.01695872 *x +  -0.001627748 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.270821389  +  0.02280435 *x +  -0.001915609 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.208666814  +  0.013038006 *x +  -0.001426277 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.19018699  +  0.011453367 *x +  2.42e-05 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.305151708  +  0.05430258 *x +  -0.001866135 *x^2}, size=1, xlim=c(0, 6 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.316049381  +  0.047485995 *x +  -0.001655464 *x^2}, size=1, xlim=c(0, 6 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.305648989  +  0.047808901 *x +  -0.001617279 *x^2}, size=1, xlim=c(0, 6 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.340442699  +  0.053605219 *x +  -0.001918748 *x^2}, size=1, xlim=c(0, 6 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.171317739  +  0.016972382 *x +  -0.000176419 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.118306868  +  0.033496233 *x +  -0.000576232 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.18628311  +  0.065670559 *x +  -0.001807713 *x^2}, size=1, xlim=c(0, 11 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.124018121  +  0.030336841 *x +  -0.000774999 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.093148945  +  0.01411747 *x +  3.03e-05 *x^2}, size=1, xlim=c(0, 7 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.122500564  +  0.033656155 *x +  -0.00076154 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.124303939  +  0.022475659 *x +  -0.000395175 *x^2}, size=1, xlim=c(0, 7 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.145399762  +  0.031374036 *x +  -0.000778897 *x^2}, size=1, xlim=c(0, 6 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.143675753  +  0.020362734 *x +  -0.001355726 *x^2}, size=1, xlim=c(0, 6 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.147866982  +  0.00952345 *x +  -0.0012212 *x^2}, size=1, xlim=c(0, 6 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.141501496  +  0.00538861 *x +  -0.001018939 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.158699733  +  0.011054617 *x +  -0.001206203 *x^2}, size=1, xlim=c(0, 6 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.166254255  +  0.016904781 *x +  -0.001413259 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.235672941  +  0.010639024 *x +  -0.00131296 *x^2}, size=1, xlim=c(0, 11 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.066290575  +  0.012227899 *x +  -0.001065694 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.199322613  +  0.013967313 *x +  -0.001127102 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.224405849  +  0.016524498 *x +  -0.001285604 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.221920595  +  0.010422075 *x +  -0.000595084 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.186889707  +  0.018288768 *x +  -0.000946616 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.189141799  +  0.015973014 *x +  -0.000896071 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.235925243  +  0.018755179 *x +  -0.001116732 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.215343212  +  0.02148559 *x +  -0.001145399 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.245393044  +  0.019511065 *x +  -0.001142569 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.223195887  +  0.031605923 *x +  -0.00144141 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.271063994  +  0.009853774 *x +  -0.001837775 *x^2}, size=1, xlim=c(0, 6 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.304939735  +  0.009329587 *x +  -0.00182886 *x^2}, size=1, xlim=c(0, 4 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.197506308  +  0.023498734 *x +  -0.002490783 *x^2}, size=1, xlim=c(0, 4 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.210232432  +  0.015816742 *x +  -0.002536284 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.077388133  +  0.007560247 *x +  -0.001474587 *x^2}, size=1, xlim=c(0, 10 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.087784757  +  0.010177828 *x +  -0.001714932 *x^2}, size=1, xlim=c(0, 10 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.158780715  +  0.000138247 *x +  0.000415305 *x^2}, size=1, xlim=c(0, 27 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.282946615  +  0.008871207 *x +  -0.000840271 *x^2}, size=1, xlim=c(0, 27 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.03341472  +  0.037417212 *x +  -0.002066182 *x^2}, size=1, xlim=c(0, 27 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.22262187  +  0.007107566 *x +  -0.001475355 *x^2}, size=1, xlim=c(0, 18 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.130798639  +  0.039854914 *x +  -0.000761886 *x^2}, size=1, xlim=c(0, 7 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.129730224  +  0.041991063 *x +  -0.000812829 *x^2}, size=1, xlim=c(0, 7 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.133538577  +  0.0286637 *x +  -0.000314958 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.133221285  +  0.036650607 *x +  -0.000616581 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.206606629  +  0.014993791 *x +  -0.001190588 *x^2}, size=1, xlim=c(0, 7 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.169838972  +  0.025115731 *x +  -0.001476843 *x^2}, size=1, xlim=c(0, 7 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.185743373  +  0.037495216 *x +  -0.001913474 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.165812088  +  0.024629841 *x +  -0.001475891 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.322452074  +  0.037698096 *x +  -0.000725417 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.220553263  +  0.023792208 *x +  7.84e-05 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.28508591  +  0.046494139 *x +  -0.000930989 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.15727932  +  0.011402923 *x +  1.12e-05 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.215567698  +  0.02113793 *x +  -0.000511491 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.282896682  +  0.018634047 *x +  -0.000471229 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.234500067  +  0.031099532 *x +  -0.000939184 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.25887029  +  0.015193302 *x +  -0.001467677 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.314332152  +  -0.001159635 *x +  -0.00096307 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.300114833  +  0.007972803 *x +  -0.001850075 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.440077859  +  0.008953618 *x +  -0.002206003 *x^2}, size=1, xlim=c(0, 5 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.401007412  +  0.006314294 *x +  -0.002004315 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.39534907  +  0.002591441 *x +  -0.001860887 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.348203078  +  0.028710656 *x +  -0.002300434 *x^2}, size=1, xlim=c(0, 5 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.327948394  +  0.006889469 *x +  -0.00191133 *x^2}, size=1, xlim=c(0, 5 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.322579102  +  0.033906501 *x +  -0.002448744 *x^2}, size=1, xlim=c(0, 5 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.162807428  +  0.01033929 *x +  -0.000990086 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.255856784  +  0.015100504 *x +  -0.00142521 *x^2}, size=1, xlim=c(0, 3 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.208078846  +  0.019881012 *x +  -0.00144848 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.256094115  +  0.020393187 *x +  -0.001511413 *x^2}, size=1, xlim=c(0, 8 ), colour='#1400E544') +
  stat_function(fun=function(x){ 0.322484612  +  0.007737378 *x +  -0.001270659 *x^2}, size=1, xlim=c(0, 8 ), colour='#4A06AC44') +
  stat_function(fun=function(x){ 0.247561942  +  0.024541204 *x +  -0.001652742 *x^2}, size=1, xlim=c(0, 3 ), colour='#800C7444') +
  stat_function(fun=function(x){ 0.241731854  +  0.031946777 *x +  -0.00197068 *x^2}, size=1, xlim=c(0, 3 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.273621638  +  0.032687607 *x +  -0.002125387 *x^2}, size=1, xlim=c(0, 3 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.298236346  +  0.028155183 *x +  -0.001922593 *x^2}, size=1, xlim=c(0, 3 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.28492165  +  0.029392281 *x +  -0.002033441 *x^2}, size=1, xlim=c(0, 3 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.420676102  +  0.031447617 *x +  -0.002369625 *x^2}, size=1, xlim=c(0, 3 ), colour='#B6123C44') +
  stat_function(fun=function(x){ 0.297726884  +  0.023773068 *x +  -0.001845754 *x^2}, size=1, xlim=c(0, 3 ), colour='#EC180444') +
  stat_function(fun=function(x){ 0.276706854  +  0.011850926 *x +  -0.001318645 *x^2}, size=1, xlim=c(0, 3 ), colour='#1400E544')
  

print(p)

print('test')
p <- ggplot(data=data.frame(x=c(0,0))) +
  ylim(c(0,1)) + xlim(c(0,28))

p <- p + stat_function(aes(y=0), fun=function(x) mean[1,'Quads']*x^2 + mean[1,'Slopes']*x + mean[1,'Intercepts'], xlim=c(0, mean[1,'experiment_length']), colour='#4A06AC') #1

p <- p + stat_function(aes(y=0), fun=function(x) -0.003370371*x^2 + 0.054212157*x + 0.53032008, xlim=c(0, 10), colour='#1400E5') #15

print(p)


# p <- ggplot(data=data.frame(x=c(0,0))) +
#   scale_y_continuous(limits=c(0,1)) +
#   scale_x_continuous(limits=c(0,28))
# 
# for (i in 1:length(mean$Intercepts))
#   p <- p + layer(stat="function", fun=function(x) mean[i, 'Quads'] * x^2 + mean[i, 'Slopes'] * x + mean[i, 'Intercepts'], xlim=c(0, mean[i, 'experiment_length']),mapping=aes(color=))
# print(p)


funcs <- list(log,function(x) x,function(x) x*log(x),function(x) x^2,  exp)
cols <-heat.colors(5,1)
p <-ggplot()+xlim(c(1,10))+ylim(c(1,10))
for(i in 1:length(funcs))
  p <- p + stat_function(aes(y=0),fun = funcs[[i]], colour=cols[i])
print(p)



ggplot(data=data.frame(x=c(0,0)), aes(x)) +
  stat_function(fun=function(x){ 0.276706854  +  0.011850926 *x +  -0.001318645 *x^2}, xlim=c(0, 3 ), colour="#1400E5" ) +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,28))



stat_function(fun=function(x) -0.00136*x^2 + 0.021503*x + 0.230929, geom='line', xlim=c(0,2), color='red')




# make plot object    
p <- qplot(
  data = data.frame(x = x, y = y), x, y, xlab = x_label, ylab = y_label, 
  enter code herexlim = x_range, main = chart_title  )

# make empty function
eq_dummy = function(x){ 0 }
d = stat_function(fun = eq_dummy)

for(i in 1 : 5){                                            
  
  # Specify Variables
  intercept = mean2[i,'Intercepts']
  slope = mean2[i,'Slopes']
  quad = mean2[i, 'Quads']
  
  # Define Curve    
  eq <- function(x) { quad*x^2 + slope*x + intercept }
  
  # Make plot object            
  composite <- stat_function(fun=eq)        
  composite = composite + d       
  
}

print(p + composite)  





#example for stack overflow

intercept <- c(0.23, 0.53, 0.41)
linear <- c(0.02, 0.05, 0.04)
quad <- c(-0.01, 0.01, 0.01)
limit <- c(5, 18, 27)
color <- c('#1400E5', '#800C74', '#EC1804')
dat <- data.frame(intercept, linear, quad, limit, color)

#solution 1
p <- ggplot(data=data.frame(x=c(0,0))) +
  ylim(c(0,1)) +
  xlim(c(0,28))

for (i in 1:length(dat$quad))
  p <- p + stat_function(aes(y=0), fun=function(x) dat[i,'quad']*x^2 + dat[i,'linear']*x + dat[i,'intercept'], xlim=c(0, dat[i,'limit']), colour=dat[i,'color'])

print(p)

#solution 2
lines <- alply(dat, 1, function(row) { 
  stat_function(fun=function(x) row['quad'] * x^2 + row['linear'] * x + row['intercept'], xlim=c(0,row['limit']), color=row['color'])  
})

ggplot(data=data.frame(x=c(0,0))) +
  lines +
  scale_y_continuous(limits=c(0,1)) +
  scale_x_continuous(limits=c(0,28))






draw <- function(row) {
  v <- stat_function(fun=function(x) {
    print(row)
    row['quad'] * x^2 + row['linear'] * x + row['intercept']
  }, xlim=c(0,row['limit']), color=row['color'])
  
  print(v)
  
  return(v)
}

p <- ggplot(data=data.frame(x=c(0,0))) +
  ylim(c(0,1)) +
  xlim(c(0,28))

for (i in 1:nrow(dat))
  p <- p + draw(dat[i,])

print(p)

