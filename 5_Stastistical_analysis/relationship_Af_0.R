
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


## Linear model with 0 and NA
###------------------------------

length(en_Af_mask_0[en_Af_mask_0==0]) # 0 cells = 0
length(en_Af_mask_0[en_Af_mask_0>0]) # 485 cells are not 0 values ~ 

lm_Af_0 <- lm(values(en_crop_0)~values(HS_Af_crop),na.action=na.omit)
summary(lm_Af_0)
#Multiple R-squared:  0.07471,  Adjusted R-squared:  0.07397 
#F-statistic: 100.3 on 1 and 1242 DF,  p-value: < 2.2e-16 ***

plot(lm_Af_0)

plot(values(en_crop_0)~values(HS_Af_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the African realm')

abline(0.022389,-0.028916 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0019287,-0.0024501,col='orange') # when 0 kept = regression in all forested area
points(x=en_Af_mask_0[en_Af_mask_0==0],y=en_Af_mask_0[en_Af_mask_0==0],col='orange')




## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_Af_0
##

## residual for lm_Af_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_Af_0$residuals))]<- lm_Af_0$residuals # replace all data-cells with res value
resid_lm_0 <- HS_Af_crop
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
# Moran I statistic standard deviate = 34.936, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#  Moran I statistic       Expectation          Variance 
#0.7534159542     -0.0008071025      0.0004660750 


# correlogram to set d2 value for the dnearneigh to calculate SARerr
##

system.time(co_lm_0 <- correlog(x,y,z_lm_0,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_0, xlim=c(100000,5000000))
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_0$correlation~co_lm_0$mean.of.class,col="red",lwd=2)
points(x=co_lm_0$x.intercept,y=0,pch=19,col="red")

# results

co_lm_0$x.intercept # 903348.7
co_lm_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_Af_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
en_0 = getValues(en_crop_0)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=903348.7,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_Af_crop)


# spatial error SAR

sar_i <- errorsarlm(en_0~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.03594866 -0.00350799 -0.00173676  0.00055804  0.07795170 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0129365  0.0073573  1.7583   0.07869
# hs          -0.0063608  0.0010911 -5.8299 5.545e-09
# 
# Lambda: 0.9657, LR test value: 435.65, p-value: < 2.22e-16
# Asymptotic standard error: 0.020089
# z-value: 48.07, p-value: < 2.22e-16
# Wald statistic: 2310.8, p-value: < 2.22e-16
# 
# Log likelihood: 4097.576 for error model
# ML residual variance (sigma squared): 7.9194e-05, (sigma: 0.0088991)
# Number of observations: 1244 
# Number of parameters estimated: 4 
# AIC: -8187.2, (AIC for lm: -7753.5)

sar_cell <- errorsarlm(en_0~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.03528854 -0.00137056 -0.00070342  0.00017747  0.05088417 
# 
# Type: error 
# Regions with no neighbours included:
#   7525 8858 9305 11833 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0067170  0.0011134  6.0326 1.613e-09
# hs          -0.0033375  0.0012045 -2.7709  0.005591
# 
# Lambda: 0.84846, LR test value: 1284, p-value: < 2.22e-16
# Asymptotic standard error: 0.014336
# z-value: 59.184, p-value: < 2.22e-16
# Wald statistic: 3502.7, p-value: < 2.22e-16
# 
# Log likelihood: 4521.757 for error model
# ML residual variance (sigma squared): 3.0773e-05, (sigma: 0.0055474)
# Number of observations: 1244 
# Number of parameters estimated: 4 
# AIC: -9035.5, (AIC for lm: -7753.5)


# Now compare how much Variation can be explained
summary(lm_Af_0)$adj.r.squared # The r_squared of the normal regression: 0.07396675

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.3480909
summary(sar_cell,Nagelkerke=T)$NK # 0.6703795 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_Af_0) #4097.576                  3879.752, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_Af_0) #  4521.757                   3879.752, p-value < 2.2e-16



### plot correlogram adding sar_Af_0 model

## residual for sar_Af_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_Af_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_Af_i

rval <- getValues(en_crop_0) # Create new raster
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

plot(co_lm_0, xlim=c(100000,5000000))
abline(h=0,lty="dotted")

#plot lm model
lines(co_lm_0$correlation~co_lm_0$mean.of.class,col="red",lwd=2)
points(x=co_lm_0$x.intercept,y=0,pch=19,col="red")

# plot sar model
lines(co_sar_i$correlation~co_sar_i$mean.of.class,col="green",lwd=2)
points(x=co_sar_i$x.intercept,y=0,pch=19,col="green")

lines(co_sar_cell$correlation~co_sar_cell$mean.of.class,col="blue",lwd=2)
points(x=co_sar_cell$x.intercept,y=0,pch=19,col="blue")
