
## IP Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: IPstraly
IP <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/IP_redefine_forest.tif') # to define IP realm
plot(IP)

IP_100km <- projectRaster(IP, HS_realm2)
IP_100km[IP_100km>=0]<-0

HS_IP_mask <- mask(HS_realm2, IP_100km)
en_IP_mask <- mask(en, IP_100km)

# put 0 in forest without small range birds
en_IP_mask_0 <- en_IP_mask
en_IP_mask_0[is.na(en_IP_mask) & !is.na(HS_IP_mask)] <- 0

# crop for faster computation
##

e <- c(5000000,180000000,-5000000,8000000)

HS_IP_crop <- crop(HS_IP_mask,e)
plot(HS_IP_crop)

en_crop_0 <- crop(en_IP_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_IP_mask,e)
plot(en_crop_NA)


## Linear model with 0 and NA
###------------------------------

length(en_IP_mask_0[en_IP_mask_0==0]) # 1296 cells = 0
length(en_IP_mask_0[en_IP_mask_0>0]) # 158 cells are not 0 values ~ 

lm_IP_NA <- lm(values(en_crop_NA)~values(HS_IP_crop),na.action=na.omit)
summary(lm_IP_NA)
# Multiple R-squared:  0.09361,  Adjusted R-squared:  0.0878 
# F-statistic: 16.11 on 1 and 156 DF,  p-value: 9.248e-05 ***
plot(values(en_crop_NA)~values(HS_IP_crop))

lm_IP_0 <- lm(values(en_crop_0)~values(HS_IP_crop),na.action=na.omit)
summary(lm_IP_0)
# Multiple R-squared:  0.09605,  Adjusted R-squared:  0.09543 
# F-statistic: 154.3 on 1 and 1452 DF,  p-value: < 2.2e-16 ***

plot(lm_IP_0)

plot(values(en_crop_0)~values(HS_IP_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the North American realm')

abline(0.0116308,-0.0051194 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0031030,-0.0038026,col='orange') # when 0 kept = regression in all forested area
points(x=HS_IP_mask[en_IP_mask_0==0],y=en_IP_mask_0[en_IP_mask_0==0],col='orange')


## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_IP_0
##

## residual for lm_IP_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_IP_0$residuals))]<- lm_IP_0$residuals # replace all data-cells with res value
resid_lm_0 <- HS_IP_crop
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
#Moran I statistic standard deviate = 20.284, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#Moran I statistic       Expectation          Variance
#     0.5198865905     -0.0010111223      0.0006594874




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

co_lm_0$x.intercept # 3436703
co_lm_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_IP_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
en_0 = getValues(en_crop_0)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=3436703,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_IP_crop)


# spatial error SAR

sar_i <- errorsarlm(en_0~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_i, na.action = na.omit,
    # zero.policy = T)

# Residuals:
       # Min         1Q     Median         3Q        Max
# -0.0941800 -0.0165720 -0.0084164 -0.0012479  0.4715132

# Type: error
# Coefficients: (asymptotic standard errors)
              # Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.3185118  0.1450559  2.1958   0.02811
# hs          -0.0233739  0.0057312 -4.0784 4.535e-05

# Lambda: 0.99022, LR test value: 289.21, p-value: < 2.22e-16
# Asymptotic standard error: 0.0067855
    # z-value: 145.93, p-value: < 2.22e-16
# Wald statistic: 21296, p-value: < 2.22e-16

# Log likelihood: 1707.282 for error model
# ML residual variance (sigma squared): 0.0020557, (sigma: 0.04534)
# Number of observations: 1022
# Number of parameters estimated: 4
# AIC: -3406.6, (AIC for lm: -3119.4)

sar_cell <- errorsarlm(en_0~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_cell, na.action = na.omit,
    # zero.policy = T)

# Residuals:
       # Min         1Q     Median         3Q        Max
# -0.1203089 -0.0124214 -0.0095670 -0.0064868  0.4655173

# Type: error
# Regions with no neighbours included:
 # 5920 6384 6914 7843 8131 8262 8393 8524 8778 9338 9431 9926 10188 10215 10344 10481 10579 10602 10757 10848 10869 10997 11001 11106 11394 11415 11634 11651 11772 11904 12855 13125 13253 13889 13891
# Coefficients: (asymptotic standard errors)
             # Estimate Std. Error z value Pr(>|z|)
# (Intercept)  0.120314   0.005733 20.9862  < 2e-16
# hs          -0.021307   0.007134 -2.9867  0.00282

# Lambda: 0.91054, LR test value: 596.49, p-value: < 2.22e-16
# Asymptotic standard error: 0.0090832
    # z-value: 100.24, p-value: < 2.22e-16
# Wald statistic: 10049, p-value: < 2.22e-16

# Log likelihood: 1860.919 for error model
# ML residual variance (sigma squared): 0.0010243, (sigma: 0.032005)
# Number of observations: 1022
# Number of parameters estimated: 4
# AIC: -3713.8, (AIC for lm: -3119.4)


# Now compare how much Variation can be explained
summary(lm_IP_0)$adj.r.squared # The r_squared of the normal regression:  0.00559088

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.2514166
summary(sar_cell,Nagelkerke=T)$NK # 0.4458014 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_IP_0) # 1707.282                  1562.676, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_IP_0) #  1860.919                   1562.676, p-value < 2.2e-16



### plot correlogram adding sar_IP_0 model

## residual for sar_IP_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_IP_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_IP_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_IP_crop
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
