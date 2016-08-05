
## SAm Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: SAmstraly
SAm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/SAm_05_realm.tif') # to define SAm realm
plot(SAm)

SAm_100km <- projectRaster(SAm, HS_realm2)
SAm_100km[SAm_100km>=0]<-0

HS_SAm_mask <- mask(HS_realm2, SAm_100km)
en_SAm_mask <- mask(en, SAm_100km)

# put 0 in forest without small range birds
en_SAm_mask_0 <- en_SAm_mask
en_SAm_mask_0[is.na(en_SAm_mask) & !is.na(HS_SAm_mask)] <- 0

# crop for faster computation
##

e <- c(-15000000,0,-15000000,5000000)

HS_SAm_crop <- crop(HS_SAm_mask,e)
plot(HS_SAm_crop)

en_crop_0 <- crop(en_SAm_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_SAm_mask,e)
plot(en_crop_NA)


## Linear model with 0 and NA
###------------------------------

length(en_crop_0[en_crop_0==0]) # 753 cells = 0
length(en_crop_0[en_crop_0>0]) # 821 cells are not 0 values ~ 

lm_SAm_NA <- lm(values(en_crop_NA)~values(HS_SAm_crop),na.action=na.omit)
summary(lm_SAm_NA)
# Multiple R-squared:  0.07774,  Adjusted R-squared:  0.07661 
# F-statistic: 69.03 on 1 and 819 DF,  p-value: 3.996e-16 ***
plot(values(en_crop_NA)~values(HS_SAm_crop))

lm_SAm_0 <- lm(values(en_crop_0)~values(HS_SAm_crop),na.action=na.omit)
summary(lm_SAm_0)
# Multiple R-squared:  0.07034,  Adjusted R-squared:  0.06975 
# F-statistic: 118.9 on 1 and 1572 DF,  p-value: < 2.2e-16

plot(lm_SAm_0)

plot(values(en_crop_0)~values(HS_SAm_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the SAmrican realm')

abline(0.022389,-0.028916 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0019287,-0.0024501,col='orange') # when 0 kept = regression in all forested area
points(x=HS_SAm_mask[en_SAm_mask_0==0],y=en_SAm_mask_0[en_SAm_mask_0==0],col='orange')



## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###-------------------------------------------------------------------

library(ncf) # For the Correlogram


# residual raster: lm_SAm_NA
##----------------------------------

## residual for lm_SAm_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_SAm_NA$residuals))]<- lm_SAm_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_SAm_crop
values(resid_lm_NA) <- rval
rm(rval) #replace our values in this new raster
names(resid_lm_NA) <- "Residuals"

plot(resid_lm_NA)

# calculate Moran's I
##---------------------------------

library(spdep) #neighbors

x = xFromCell(resid_lm_NA,1:ncell(resid_lm_NA)) # take x coordinates
y = yFromCell(resid_lm_NA,1:ncell(resid_lm_NA)) # take y coordinates
z_lm_NA = getValues(resid_lm_NA) # and the values of course

##create a list of neighbors
nb <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F) # max dist = cell size for first setting
w_mtx <-nb2listw(nb,style="W",zero.policy=T) # weight

##Moran's I
#we need data + weights

moran.test(values(en_crop_NA),w_mtx, na.action=na.omit, zero.policy=T)
#Moran I statistic standard deviate = 29.045, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#  Moran I statistic       Expectation          Variance 
#0.8987384359     -0.0012562814      0.0009601129 


# correlogram to set d2 value for the dnearneigh to calculate SARerr
##

system.time(co_lm_NA <- correlog(x,y,z_lm_NA,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_NA)
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_NA$correlation~co_lm_NA$mean.of.class,col="red",lwd=2)
points(x=co_lm_NA$x.intercept,y=0,pch=19,col="red")

# results

co_lm_NA$x.intercept # 3526477
co_lm_NA$correlation 


## SAR model
###------------------------------------------------------------###

# create a list of neighbors
##-----------------------------------

library(spdep)

x = xFromCell(en_crop_NA,1:ncell(en_crop_NA))
y = yFromCell(en_crop_NA,1:ncell(en_crop_NA))
en_NA = getValues(en_crop_NA)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=3526477,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_SAm_crop)


# spatial error SAR
##----------------------------------

sar_i <- errorsarlm(en_NA~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0527913 -0.0158697 -0.0057684  0.0035497  0.3046901 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error  z value Pr(>|z|)
# (Intercept)  0.0255410  0.0372264   0.6861   0.4926
# hs          -0.0445527  0.0040893 -10.8949   <2e-16
# 
# Lambda: 0.96783, LR test value: 102.95, p-value: < 2.22e-16
# Asymptotic standard error: 0.022379
# z-value: 43.247, p-value: < 2.22e-16
# Wald statistic: 1870.3, p-value: < 2.22e-16
# 
# Log likelihood: 1601.413 for error model
# ML residual variance (sigma squared): 0.0011757, (sigma: 0.034288)
# Number of observations: 821 
# Number of parameters estimated: 4 
# AIC: -3194.8, (AIC for lm: -3093.9)


sar_cell <- errorsarlm(en_NA~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0884787 -0.0042359 -0.0014087  0.0024896  0.1095287 
# 
# Type: error 
# Regions with no neighbours included:
#   2466 2643 3228 3393 3529 3691 4293 7487 8110 8714 9021 9025 9293 9306 9627 10514 10826 11130 11423 12924 15341 15644 15796 17010 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0227990  0.0027485  8.2951 < 2.2e-16
# hs          -0.0132346  0.0034325 -3.8557 0.0001154
# 
# Lambda: 0.86061, LR test value: 1282.1, p-value: < 2.22e-16
# Asymptotic standard error: 0.012115
# z-value: 71.038, p-value: < 2.22e-16
# Wald statistic: 5046.3, p-value: < 2.22e-16
# 
# Log likelihood: 2190.97 for error model
# ML residual variance (sigma squared): 0.00018851, (sigma: 0.01373)
# Number of observations: 821 
# Number of parameters estimated: 4 
# AIC: -4373.9, (AIC for lm: -3093.9)


# Now compare how much Variation can be explained
summary(lm_SAm_NA)$adj.r.squared # The r_squared of the normal regression: 0.07661189

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.1864332
summary(sar_cell,Nagelkerke=T)$NK # 0.8065081 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_SAm_NA) #  1601.413                    1549.936  , < 2.2e-16

LR.sarlm(sar_cell,lm_SAm_NA) #   2190.970                    1549.936 , p-value < 2.2e-16



### plot correlogram adding sar_SAm_NA model
##-------------------------------------------

## residual for sar_SAm_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_SAm_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_SAm_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_SAm_crop
values(resid_sar_cell) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_cell) <- "Residuals"

# to get values
x = xFromCell(resid_sar_cell,1:ncell(resid_sar_cell)) # take x coordinates
y = yFromCell(resid_sar_cell,1:ncell(resid_sar_cell)) # take y coordinates

z_sar_i = getValues(resid_sar_i) 
z_sar_cell = getValues(resid_sar_cell) 

# sar correlogram
system.time(co_sar_i <- correlog(x,y,z_sar_i,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

system.time(co_sar_cell <- correlog(x,y,z_sar_cell,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

## plot correlogram

plot(co_lm_NA, xlim=c(100000,5000000))
abline(h=0,lty="dotted")

#plot lm model
lines(co_lm_NA$correlation~co_lm_NA$mean.of.class,col="red",lwd=2)
points(x=co_lm_NA$x.intercept,y=0,pch=19,col="red")

# plot sar model
lines(co_sar_i$correlation~co_sar_i$mean.of.class,col="green",lwd=2)
points(x=co_sar_i$x.intercept,y=0,pch=19,col="green")

lines(co_sar_cell$correlation~co_sar_cell$mean.of.class,col="blue",lwd=2)
points(x=co_sar_cell$x.intercept,y=0,pch=19,col="blue")

# test pr plotter sar model

# (Intercept)  0.0227990  0.0027485  8.2951 < 2.2e-16
# hs          -0.0132346  0.0034325 -3.8557 0.0001154

# (Intercept)  0.0255410  0.0372264   0.6861   0.4926
# hs          -0.0445527  0.0040893 -10.8949   <2e-16

#plot(z_sar_cell~hs)
#abline(0.0227990,-0.0132346 , col="red") # cell
#abline(0.0255410,-0.0445527 , col="orange") # i