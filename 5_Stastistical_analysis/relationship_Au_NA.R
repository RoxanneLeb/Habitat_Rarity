

## Au Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: Australy
Au <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Au_02_realm.tif') # to define AU realm
plot(Au)

# redifine resolution : 100 km
Au_100km <- projectRaster(Au, HS_realm2)
Au_100km[Au_100km>=0]<-0

# redifine study area
HS_Au_mask <- mask(HS_realm2, Au_100km)
en_Au_mask <- mask(en, Au_100km)

# put 0 in area without small range birds
en_Au_mask_0 <- en_Au_mask
en_Au_mask_0[is.na(en_Au_mask) & !is.na(HS_Au_mask)] <- 0

# crop for faster computation

e <- c(5000000,180000000,-8000000,0)

HS_Au_crop <- crop(HS_Au_mask,e)
plot(HS_Au_crop)

en_crop_0 <- crop(en_Au_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_Au_mask,e)
plot(en_crop_NA)


## Linear model 
###------------------------------

length(en_Au_mask_0[en_Au_mask_0==0]) # 172 cells = 0
length(en_Au_mask_0[en_Au_mask_0>0]) # 173 cells are not 0 values ~ 50%

## lm_Au_NA

lm_Au_NA <- lm(values(en_crop_NA)~values(HS_Au_crop),na.action=na.omit)
summary(lm_Au_NA)
#Multiple R-squared:  0.08765,  Adjusted R-squared:  0.08232 
#F-statistic: 16.43 on 1 and 171 DF,  p-value: 7.658e-05
plot(values(en_crop_NA)~values(HS_Au_crop))

## lm_Au_0

lm_Au_0 <- lm(values(en_crop_0)~values(HS_Au_crop),na.action=na.omit)
summary(lm_Au_0)
#Multiple R-squared:  0.04705,  Adjusted R-squared:  0.04427 
#F-statistic: 16.93 on 1 and 343 DF,  p-value: 4.853e-05

plot(lm_Au_0)

plot(values(en_crop_NA)~values(HS_Au_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the Australian realm')

abline(0.020669,-0.017539,col='orange') # when 0 kept = regression in all forested area
abline(0.037467,-0.029049 ,col='red') # when 0 removed = regression just inside sites where small species occur
points(x=HS_Au_mask[en_Au_mask_0==0],y=en_Au_mask_0[en_Au_mask_0==0],col='orange')


## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_Au_NA
##

## residual for lm_Au_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_Au_NA$residuals))]<- lm_Au_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_Au_crop
values(resid_lm_NA) <- rval
rm(rval) #replace our values in this new raster
names(resid_lm_NA) <- "Residuals"

resid_lm_NA

# calculate Moran's I
##-------------------

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
# Moran I statistic standard deviate = 8.5999, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
#  Moran I statistic       Expectation          Variance 
#  0.652569057      -0.006211180       0.005868066 



# correlogram to set d2 value for the dnearneigh to calculate SARerr
system.time(co_lm_NA <- correlog(x,y,z_lm_NA,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_NA, xlim=c(100000,5000000))
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_NA$correlation~co_lm_NA$mean.of.class,col="red",lwd=2)
points(x=co_lm_NA$x.intercept,y=0,pch=19,col="red")

# results

co_lm_NA$x.intercept # 860506.7
co_lm_NA$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_Au_NA
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_NA,1:ncell(en_crop_NA))
y = yFromCell(en_crop_NA,1:ncell(en_crop_NA))
en_NA = getValues(en_crop_NA)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=860506.7,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_Au_crop)


# spatial error SAR

sar_i <- errorsarlm(en_NA~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_i, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0378588 -0.0095764 -0.0052250  0.0022598  0.1004033 
# 
# Type: error 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.0304584  0.0092649  3.2875 0.001011
# hs          -0.0093246  0.0074833 -1.2460 0.212746
# 
# Lambda: 0.82181, LR test value: 78.14, p-value: < 2.22e-16
# Asymptotic standard error: 0.062569
# z-value: 13.134, p-value: < 2.22e-16
# Wald statistic: 172.51, p-value: < 2.22e-16
# 
# Log likelihood: 426.1997 for error model
# ML residual variance (sigma squared): 0.0004006, (sigma: 0.020015)
# Number of observations: 173 
# Number of parameters estimated: 4 
# AIC: -844.4, (AIC for lm: -768.26)

sar_cell <- errorsarlm(en_NA~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_cell, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min         1Q     Median         3Q        Max 
# -0.0255295 -0.0081346 -0.0044333  0.0032389  0.1128904 
# 
# Type: error 
# Regions with no neighbours included:
#   1899 2162 2414 2543 3982 5142 5144 5401 5693 6722 7233 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0319980  0.0047534  6.7316 1.678e-11
# hs          -0.0201913  0.0080571 -2.5060   0.01221
# 
# Lambda: 0.59501, LR test value: 81.928, p-value: < 2.22e-16
# Asymptotic standard error: 0.051695
# z-value: 11.51, p-value: < 2.22e-16
# Wald statistic: 132.48, p-value: < 2.22e-16
# 
# Log likelihood: 428.0932 for error model
# ML residual variance (sigma squared): 0.00034722, (sigma: 0.018634)
# Number of observations: 173 
# Number of parameters estimated: 4 
# AIC: -848.19, (AIC for lm: -768.26)


# Now compare how much Variation can be explained
summary(lm_Au_NA)$adj.r.squared # The r_squared of the normal regression: 0.08231617

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.4192365
summary(sar_cell,Nagelkerke=T)$NK # 0.4318119 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_Au_NA) #sar_e: 426.1997 lm_Au_NA: 387.1295, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_Au_NA) #  428.0932                   387.1295 , p-value < 2.2e-16



### plot correlogram adding sar_Au_NA model

## residual for sar_Au_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_Au_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_Au_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_Au_crop
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
