
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
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_NAm_0
##

## residual for lm_NAm_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_NAm_0$residuals))]<- lm_NAm_0$residuals # replace all data-cells with res value
resid_lm_0 <- HS_NAm_crop
values(resid_lm_0) <- rval
rm(rval) #replace our values in this new raster
names(resid_lm_0) <- "Residuals"

plot(resid_lm_0)

# calculate Moran's I
##-------------------

library(spdep) #neighbors

x = xFromCell(resid_lm_0,1:ncell(resid_lm_0)) # take x coordinates
y = yFromCell(resid_lm_0,1:ncell(resid_lm_0)) # take y coordinates
z_lm_0 = getValues(resid_lm_0) # and the values of course

##create a list of neighbors
nb <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F) # max dist = cell size for first setting
w_mtx <-nb2listw(nb,style="W",zero.policy=T) # weight

##Moran's I
#we need data + weights

moran.test(values(en_crop_0),w_mtx, na.action=na.omit, zero.policy=T)
# Moran I statistic standard deviate = 29.889, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
#  Moran I statistic       Expectation          Variance 
# 0.6028937849     -0.0006930007      0.0004078007 



# correlogram to set d2 value for the dnearneigh to calculate SARerr
##

system.time(co_lm_0 <- correlog(x,y,z_lm_0,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_0)
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_0$correlation~co_lm_0$mean.of.class,col="red",lwd=2)
points(x=co_lm_0$x.intercept,y=0,pch=19,col="red")

# results

co_lm_0$x.intercept # 843937.6
co_lm_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_NAm_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
en_0 = getValues(en_crop_0)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=843937.6,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_NAm_crop)


# spatial error SAR

sar_i <- errorsarlm(en_0~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.00798020 -0.00105072 -0.00047116 -0.00011835  0.02693608 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate  Std. Error z value  Pr(>|z|)
# (Intercept)  0.00428201  0.00191154  2.2401 0.0250851
# hs          -0.00149477  0.00042204 -3.5417 0.0003975
# 
# Lambda: 0.95828, LR test value: 255.62, p-value: < 2.22e-16
# Asymptotic standard error: 0.02148
# z-value: 44.613, p-value: < 2.22e-16
# Wald statistic: 1990.4, p-value: < 2.22e-16
# 
# Log likelihood: 6353.533 for error model
# ML residual variance (sigma squared): 9.2155e-06, (sigma: 0.0030357)
# Number of observations: 1454 
# Number of parameters estimated: 4 
# AIC: -12699, (AIC for lm: -12445)

sar_cell <- errorsarlm(en_0~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.01069375 -0.00052977 -0.00032663 -0.00018111  0.02796832 
# 
# Type: error 
# Regions with no neighbours included:
#   2457 2624 2924 3698 4029 4035 4179 4599 4910 7898 
# Coefficients: (asymptotic standard errors) 
# Estimate  Std. Error z value  Pr(>|z|)
# (Intercept)  0.00210545  0.00030063  7.0035 2.497e-12
# hs          -0.00150342  0.00044642 -3.3677  0.000758
# 
# Lambda: 0.71351, LR test value: 747.47, p-value: < 2.22e-16
# Asymptotic standard error: 0.019837
# z-value: 35.969, p-value: < 2.22e-16
# Wald statistic: 1293.8, p-value: < 2.22e-16
# 
# Log likelihood: 6599.458 for error model
# ML residual variance (sigma squared): 5.5983e-06, (sigma: 0.0023661)
# Number of observations: 1454 
# Number of parameters estimated: 4 
# AIC: -13191, (AIC for lm: -12445)


# Now compare how much Variation can be explained
summary(lm_NAm_0)$adj.r.squared # The r_squared of the normal regression:  0.09542534

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.2417804
summary(sar_cell,Nagelkerke=T)$NK # 0.4593894 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_NAm_0) #6353.533                   6225.724, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_NAm_0) #  6599.458                   6225.724, p-value < 2.2e-16



### plot correlogram adding sar_NAm_0 model

## residual for sar_NAm_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_NAm_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_NAm_i

rval <- getValues(en_crop_0) # Create new raster
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

plot(co_lm_0)
abline(h=0,lty="dotted")

#plot lm model
lines(co_lm_0$correlation~co_lm_0$mean.of.class,col="red",lwd=2)
points(x=co_lm_0$x.intercept,y=0,pch=19,col="red")

# plot sar model
lines(co_sar_i$correlation~co_sar_i$mean.of.class,col="green",lwd=2)
points(x=co_sar_i$x.intercept,y=0,pch=19,col="green")

lines(co_sar_cell$correlation~co_sar_cell$mean.of.class,col="blue",lwd=2)
points(x=co_sar_cell$x.intercept,y=0,pch=19,col="blue")
