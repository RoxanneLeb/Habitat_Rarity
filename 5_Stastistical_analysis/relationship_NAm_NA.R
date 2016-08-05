
## NAm Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: NAmstraly
NAm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/NAm_06_realm.tif') # to define NAm realm
plot(NAm)

NAm_100km <- projectRaster(NAm, HS_realm2)
NAm_100km[NAm_100km>=0]<-0

HS_NAm_mask <- mask(HS_realm2, NAm_100km)
en_NAm_mask <- mask(en, NAm_100km)

# put 0 in forest without small range birds
en_NAm_mask_0 <- en_NAm_mask
en_NAm_mask_0[is.na(en_NAm_mask) & !is.na(HS_NAm_mask)] <- 0

# crop for faster computation
##

e <- c(-15000000,0,0,15000000)

HS_NAm_crop <- crop(HS_NAm_mask,e)
plot(HS_NAm_crop)

en_crop_0 <- crop(en_NAm_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_NAm_mask,e)
plot(en_crop_NA)


## Linear model with 0 and NA
###------------------------------

length(en_NAm_mask_0[en_NAm_mask_0==0]) # 1296 cells = 0
length(en_NAm_mask_0[en_NAm_mask_0>0]) # 158 cells are not 0 values ~ 

lm_NAm_NA <- lm(values(en_crop_NA)~values(HS_NAm_crop),na.action=na.omit)
summary(lm_NAm_NA)
# Multiple R-squared:  0.09361,  Adjusted R-squared:  0.0878 
# F-statistic: 16.11 on 1 and 156 DF,  p-value: 9.248e-05 ***
plot(values(en_crop_NA)~values(HS_NAm_crop))

lm_NAm_0 <- lm(values(en_crop_0)~values(HS_NAm_crop),na.action=na.omit)
summary(lm_NAm_0)
# Multiple R-squared:  0.09605,  Adjusted R-squared:  0.09543 
# F-statistic: 154.3 on 1 and 1452 DF,  p-value: < 2.2e-16 ***

plot(lm_NAm_0)

plot(values(en_crop_0)~values(HS_NAm_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the North American realm')

abline(0.0116308,-0.0051194 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0031030,-0.0038026,col='orange') # when 0 kept = regression in all forested area
points(x=HS_NAm_mask[en_NAm_mask_0==0],y=en_NAm_mask_0[en_NAm_mask_0==0],col='orange')



## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###-------------------------------------------------------------------

library(ncf) # For the Correlogram


# residual raster: lm_NAm_NA
##----------------------------------

## residual for lm_NAm_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_NAm_NA$residuals))]<- lm_NAm_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_NAm_crop
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
# Moran I statistic standard deviate = 6.706, p-value = 1e-11
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.560738421      -0.006944444       0.007166224 


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

co_lm_NA$x.intercept # 859754.2
co_lm_NA$correlation 


## SAR model
###------------------------------------------------------------###

# create a list of neighbors
##-----------------------------------

library(spdep)

x = xFromCell(en_crop_NA,1:ncell(en_crop_NA))
y = yFromCell(en_crop_NA,1:ncell(en_crop_NA))
en_NA = getValues(en_crop_NA)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=859754.2,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_NAm_crop)


# spatial error SAR
##----------------------------------

sar_i <- errorsarlm(en_NA~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.00748874 -0.00238845 -0.00100104  0.00072496  0.01835387 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate  Std. Error z value  Pr(>|z|)
# (Intercept)  0.01191898  0.00081035 14.7085 < 2.2e-16
# hs          -0.00580640  0.00169192 -3.4318 0.0005995
# 
# Lambda: 0.49587, LR test value: 11.018, p-value: 0.00090244
# Asymptotic standard error: 0.13531
# z-value: 3.6648, p-value: 0.00024755
# Wald statistic: 13.431, p-value: 0.00024755
# 
# Log likelihood: 639.2828 for error model
# ML residual variance (sigma squared): 1.7569e-05, (sigma: 0.0041916)
# Number of observations: 158 
# Number of parameters estimated: 4 
# AIC: -1270.6, (AIC for lm: -1261.5)


sar_cell <- errorsarlm(en_NA~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.00763335 -0.00188566 -0.00095883  0.00043239  0.01676688 
# 
# Type: error 
# Regions with no neighbours included:
#   1719 2324 2924 2926 3530 3540 3689 5071 5219 5221 5671 6690 7898 
# Coefficients: (asymptotic standard errors) 
# Estimate  Std. Error z value  Pr(>|z|)
# (Intercept)  0.01180481  0.00067087 17.5961 < 2.2e-16
# hs          -0.00481337  0.00160420 -3.0005  0.002695
# 
# Lambda: 0.49558, LR test value: 48.962, p-value: 2.6102e-12
# Asymptotic standard error: 0.059808
# z-value: 8.2862, p-value: 2.2204e-16
# Wald statistic: 68.662, p-value: < 2.22e-16
# 
# Log likelihood: 658.2547 for error model
# ML residual variance (sigma squared): 1.242e-05, (sigma: 0.0035242)
# Number of observations: 158 
# Number of parameters estimated: 4 
# AIC: -1308.5, (AIC for lm: -1261.5)


# Now compare how much Variation can be explained
summary(lm_NAm_NA)$adj.r.squared # The r_squared of the normal regression: 0.08780258

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.154664
summary(sar_cell,Nagelkerke=T)$NK # 0.3351354 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_NAm_NA) #  639.2828                    633.7739  , = 0.0009024

LR.sarlm(sar_cell,lm_NAm_NA) #   658.2547                    633.7739  , p-value = 2.61e-12



### plot correlogram adding sar_NAm_NA model
##-------------------------------------------

## residual for sar_NAm_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_NAm_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_NAm_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_NAm_crop
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
