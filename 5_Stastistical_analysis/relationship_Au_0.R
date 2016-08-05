

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


## Linear model with 0 
###------------------------------

length(en_Au_mask_0[en_Au_mask_0==0]) # 172 cells = 0
length(en_Au_mask_0[en_Au_mask_0>0]) # 173 cells are not 0 values ~ 50%

lm_Au_0 <- lm(values(en_crop_0)~values(HS_Au_crop),na.action=na.omit)
summary(lm_Au_0)
#Multiple R-squared:  0.04705,  Adjusted R-squared:  0.04427 
#F-statistic: 16.93 on 1 and 343 DF,  p-value: 4.853e-05

plot(lm_Au_0)

plot(values(en_crop_0)~values(HS_Au_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the Australian realm')

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


# lm_Au_0
##

## residual for lm_Au_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_Au_0$residuals))]<- lm_Au_0$residuals # replace all data-cells with res value
resid_0 <- HS_Au_crop
values(resid_0) <- rval
rm(rval) #replace our values in this new raster
names(resid_0) <- "Residuals"

resid_0
plot(resid_0)


# calculate Moran's I
##-------------------

x = xFromCell(resid_0,1:ncell(resid_0)) # take x coordinates
y = yFromCell(resid_0,1:ncell(resid_0)) # take y coordinates
z_0 = getValues(resid_0) # and the values of course

##create a list of neighbors
nb <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F) # max dist = cell size for first setting
w_mtx <-nb2listw(nb,style="W",zero.policy=T) # weight

##Moran's I
#we need data + weights

moran.test(values(en_crop_0),w_mtx, na.action=na.omit, zero.policy=T)
# Moran I statistic standard deviate = 14.087, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
#   Moran I statistic       Expectation          Variance 
# 0.684299722      -0.003048780       0.002380648 


# correlogram to set d2 value for the dnearneigh to calculate SARerr
system.time(co_0 <- correlog(x,y,z_0,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot correlogram
plot(co_0, xlim=c(100000,5000000))
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_0$correlation~co_0$mean.of.class,col="red",lwd=2)
points(x=co_0$x.intercept,y=0,pch=19,col="red")

# results

co_0$x.intercept # 3174985
co_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_Au_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
z_0 = getValues(en_crop_0)

nb <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
wei <-nb2listw(nb,style="W",zero.policy=T) 


nb_i <- dnearneigh(cbind(x,y),d1=0,d2=3174985,longlat=F)# Get neighbourlist of interactions with a distance unit d2=intercept
wei_i <-nb2listw(nb_i,style="W",zero.policy=T) 

# fit the spatial error SAR
hs <- values(HS_Au_crop)
sar_Au_0 <- errorsarlm(z_0 ~ hs, listw=wei_i, na.action=na.omit, zero.policy=T)

sar_Au_0_res <- errorsarlm(z_0 ~ hs, listw=wei, na.action=na.omit, zero.policy=T)

# We use the generated z values and weights as input. Nodata values are excluded and zeros are given to boundary errors

summary(sar_Au_0) #Lambda: 0.94331, LR test value: 29.778, p-value: 4.8434e-08

summary(sar_Au_0_res) #Lambda: 0.64887, LR test value: 199.16, p-value: < 2.22e-16
# 
# Call:errorsarlm(formula = z_0 ~ hs, listw = wei, na.action = na.omit,     zero.policy = T)
# 
# Residuals:
#   Min          1Q      Median          3Q         Max 
# -0.05278299 -0.00456608 -0.00292870  0.00065106  0.11930700 
# 
# Type: error 
# Regions with no neighbours included:
#   2414 2423 2430 2543 2686 2929 2936 3185 3187 3313 3441 3828 3982 5010 5693 7233 
# Coefficients: (asymptotic standard errors) 
# Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0145932  0.0029461  4.9535 7.291e-07
# hs          -0.0078627  0.0049780 -1.5795    0.1142
# 
# Lambda: 0.64887, LR test value: 199.16, p-value: < 2.22e-16
# Asymptotic standard error: 0.037365
# z-value: 17.366, p-value: < 2.22e-16
# Wald statistic: 301.57, p-value: < 2.22e-16
# 
# Log likelihood: 932.6927 for error model
# ML residual variance (sigma squared): 0.00021883, (sigma: 0.014793)
# Number of observations: 345 
# Number of parameters estimated: 4 
# AIC: -1857.4, (AIC for lm: -1660.2)

# Now compare how much Variation can be explained
summary(lm_Au_0)$adj.r.squared # The r_squared of the normal regression: 0.04426793

summary(sar_Au_0,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the sar_Au_0: 0.12585, with sar_Au_0_res: 0.4649825

# Finally do a likelihood ratio test
LR.sarlm(sar_Au_0,lm_Au_0) 
#sar_Au_0 :848.0034            lm: 833.1142 
#sar_Au_0_res : 932.6927       lm: 833.1142  p-value < 2.2e-16


### plot correlogram adding sar_Au_0 model

## residual for sar_Au_0

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_Au_0$residuals))]<- sar_Au_0$residuals # replace all data-cells with res value
resid_sar_0 <- HS_Au_crop
values(resid_sar_0) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_0) <- "Residuals"


## residual for sar_Au_0_res

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_Au_0_res$residuals))]<- sar_Au_0_res$residuals # replace all data-cells with res value
resid_sar_0_res <- HS_Au_crop
values(resid_sar_0_res) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_0_res) <- "Residuals"

# to get values
x = xFromCell(resid_sar_0,1:ncell(resid_sar_0)) # take x coordinates
y = yFromCell(resid_sar_0,1:ncell(resid_sar_0)) # take y coordinates

z_sar_0 = getValues(resid_sar_0) 
z_sar_0_res = getValues(resid_sar_0_res) 

# sar correlogram
system.time(co_sar_0 <- correlog(x,y,z_sar_0,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

system.time(co_sar_0_res <- correlog(x,y,z_sar_0_res,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

## plot correlogram

plot(co_0, xlim=c(100000,5000000))
abline(h=0,lty="dotted")

#plot lm model
lines(co_0$correlation~co_0$mean.of.class,col="red",lwd=2)
points(x=co_0$x.intercept,y=0,pch=19,col="red")

# plot sar model
lines(co_sar_0$correlation~co_sar_0$mean.of.class,col="green",lwd=2)
points(x=co_sar_0$x.intercept,y=0,pch=19,col="green")

lines(co_sar_0_res$correlation~co_sar_0_res$mean.of.class,col="blue",lwd=2)
points(x=co_sar_0_res$x.intercept,y=0,pch=19,col="blue")
