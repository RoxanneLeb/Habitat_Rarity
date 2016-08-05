
## Af Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: Afstraly
Af <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Af_03_realm.tif') # to define Af realm
plot(Af)

Af_100km <- projectRaster(Af, HS_realm2)
Af_100km[Af_100km>=0]<-0

HS_Af_mask <- mask(HS_realm2, Af_100km)
en_Af_mask <- mask(en, Af_100km)

# put 0 in forest without small range birds
en_Af_mask_0 <- en_Af_mask
en_Af_mask_0[is.na(en_Af_mask) & !is.na(HS_Af_mask)] <- 0

# crop for faster computation
##

e <- c(-5000000,6000000,-5000000,10000000)

HS_Af_crop <- crop(HS_Af_mask,e)
plot(HS_Af_crop)

en_crop_0 <- crop(en_Af_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_Af_mask,e)
plot(en_crop_NA)


## Linear model with 0 and NA
###------------------------------

length(en_Af_mask_0[en_Af_mask_0==0]) # 759 cells = 0
length(en_Af_mask_0[en_Af_mask_0>0]) # 485 cells are not 0 values ~ 

lm_Af_NA <- lm(values(en_crop_NA)~values(HS_Af_crop),na.action=na.omit)
summary(lm_Af_NA)
# Multiple R-squared:  0.06797,  Adjusted R-squared:  0.06604 
# F-statistic: 35.22 on 1 and 483 DF,  p-value: 5.613e-09 ***
plot(values(en_crop_NA)~values(HS_Af_crop))

lm_Af_0 <- lm(values(en_crop_0)~values(HS_Af_crop),na.action=na.omit)
summary(lm_Af_0)
# Multiple R-squared:  0.07471,  Adjusted R-squared:  0.07397 
# F-statistic: 100.3 on 1 and 1242 DF,  p-value: < 2.2e-16

plot(lm_Af_0)

plot(values(en_crop_0)~values(HS_Af_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the African realm')

abline(0.022389,-0.028916 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0019287,-0.0024501,col='orange') # when 0 kept = regression in all forested area
points(x=HS_Af_mask[en_Af_mask_0==0],y=en_Af_mask_0[en_Af_mask_0==0],col='orange')



## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###-------------------------------------------------------------------

library(ncf) # For the Correlogram


# residual raster: lm_Af_NA
##----------------------------------

## residual for lm_Af_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_Af_NA$residuals))]<- lm_Af_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_Af_crop
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
# Moran I statistic standard deviate = 17.492, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.723930120      -0.002114165       0.001722885 


# correlogram to set d2 value for the dnearneigh to calculate SARerr
##

system.time(co_lm_NA <- correlog(x,y,z_lm_NA,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_NA, xlim=c(100000,5000000))
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_NA$correlation~co_lm_NA$mean.of.class,col="red",lwd=2)
points(x=co_lm_NA$x.intercept,y=0,pch=19,col="red")

# results

co_lm_NA$x.intercept # 975060
co_lm_NA$correlation 


## SAR model
###------------------------------------------------------------###

# create a list of neighbors
##-----------------------------------

library(spdep)

x = xFromCell(en_crop_NA,1:ncell(en_crop_NA))
y = yFromCell(en_crop_NA,1:ncell(en_crop_NA))
en_NA = getValues(en_crop_NA)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=860506.7,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_Af_crop)


# spatial error SAR
##----------------------------------

sar_i <- errorsarlm(en_NA~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0306057 -0.0061735 -0.0025575  0.0016732  0.0728585 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0170754  0.0045345  3.7657 0.0001661
# hs          -0.0074520  0.0022953 -3.2467 0.0011676
# 
# Lambda: 0.87927, LR test value: 181.17, p-value: < 2.22e-16
# Asymptotic standard error: 0.043888
# z-value: 20.035, p-value: < 2.22e-16
# Wald statistic: 401.39, p-value: < 2.22e-16
# 
# Log likelihood: 1456.794 for error model
# ML residual variance (sigma squared): 0.0001391, (sigma: 0.011794)
# Number of observations: 485 
# Number of parameters estimated: 4 
# AIC: -2905.6, (AIC for lm: -2726.4)

sar_cell <- errorsarlm(en_NA~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0330886 -0.0032479 -0.0011501  0.0016589  0.0560532 
# 
# Type: error 
# Regions with no neighbours included:
#   7525 7635 7697 7856 8858 9273 9305 9385 9720 12957 13507 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0135561  0.0015627  8.6750 < 2.2e-16
# hs          -0.0086637  0.0024096 -3.5955 0.0003238
# 
# Lambda: 0.70751, LR test value: 350.28, p-value: < 2.22e-16
# Asymptotic standard error: 0.02693
# z-value: 26.272, p-value: < 2.22e-16
# Wald statistic: 690.21, p-value: < 2.22e-16
# 
# Log likelihood: 1541.352 for error model
# ML residual variance (sigma squared): 7.975e-05, (sigma: 0.0089303)
# Number of observations: 485 
# Number of parameters estimated: 4 
# AIC: -3074.7, (AIC for lm: -2726.4)


# Now compare how much Variation can be explained
summary(lm_Af_NA)$adj.r.squared # The r_squared of the normal regression: 0.06604179

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.3584884
summary(sar_cell,Nagelkerke=T)$NK # 0.5473425 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_Af_NA) # 1456.794                   1366.211 , p-value < 2.2e-16

LR.sarlm(sar_cell,lm_Af_NA) #  1541.352                   1366.211 , p-value < 2.2e-16



### plot correlogram adding sar_Af_NA model
##-------------------------------------------

## residual for sar_Af_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_Af_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_Af_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_Af_crop
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

plot(co_0, xlim=c(100000,5000000))
abline(h=0,lty="dotted")

#plot lm model
lines(co_0$correlation~co_0$mean.of.class,col="red",lwd=2)
points(x=co_0$x.intercept,y=0,pch=19,col="red")

# plot sar model
lines(co_sar_i$correlation~co_sar_i$mean.of.class,col="green",lwd=2)
points(x=co_sar_i$x.intercept,y=0,pch=19,col="green")

lines(co_sar_cell$correlation~co_sar_cell$mean.of.class,col="blue",lwd=2)
points(x=co_sar_cell$x.intercept,y=0,pch=19,col="blue")
