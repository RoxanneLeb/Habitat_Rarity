
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

length(en_IP_mask_0[en_IP_mask_0==0]) # 463 cells = 0
length(en_IP_mask_0[en_IP_mask_0>0]) # 564 cells are not 0 values ~ 

lm_IP_NA <- lm(values(en_crop_NA)~values(HS_IP_crop),na.action=na.omit)
summary(lm_IP_NA)
# Multiple R-squared:  0.01438,  Adjusted R-squared:  0.01261 
# F-statistic: 8.126 on 1 and 557 DF,  p-value: 0.004524 **
plot(values(en_crop_NA)~values(HS_IP_crop))

lm_IP_0 <- lm(values(en_crop_0)~values(HS_IP_crop),na.action=na.omit)
summary(lm_IP_0)
# Multiple R-squared:  0.006565,  Adjusted R-squared:  0.005591 
# F-statistic:  6.74 on 1 and 1020 DF,  p-value: 0.009561 **

plot(lm_IP_0)

plot(values(en_crop_0)~values(HS_IP_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the Indo-Pacific realm')

abline(0.0116308,-0.0051194 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0031030,-0.0038026,col='orange') # when 0 kept = regression in all forested area
points(x=HS_IP_mask[en_IP_mask_0==0],y=en_IP_mask_0[en_IP_mask_0==0],col='orange')



## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###-------------------------------------------------------------------

library(ncf) # For the Correlogram


# residual raster: lm_IP_NA
##----------------------------------

## residual for lm_IP_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_IP_NA$residuals))]<- lm_IP_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_IP_crop
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
#Moran I statistic standard deviate = 12.737, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#Moran I statistic       Expectation          Variance
#      0.503149817      -0.001901141       0.001572259



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

co_lm_NA$x.intercept # 4434830
co_lm_NA$correlation 


## SAR model
###------------------------------------------------------------###

# create a list of neighbors
##-----------------------------------

library(spdep)

x = xFromCell(en_crop_NA,1:ncell(en_crop_NA))
y = yFromCell(en_crop_NA,1:ncell(en_crop_NA))
en_NA = getValues(en_crop_NA)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=4434830,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_IP_crop)


# spatial error SAR
##----------------------------------

sar_i <- errorsarlm(en_NA~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_i, na.action = na.omit,
    # zero.policy = T)

# Residuals:
      # Min        1Q    Median        3Q       Max
# -0.082025 -0.026745 -0.015261 -0.000559  0.494033

# Type: error
# Coefficients: (asymptotic standard errors)
             # Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.186199   0.111816  1.6652   0.09587
# hs          -0.046474   0.010382 -4.4765 7.589e-06

# Lambda: 0.97747, LR test value: 111.77, p-value: < 2.22e-16
# Asymptotic standard error: 0.015751
    # z-value: 62.056, p-value: < 2.22e-16
# Wald statistic: 3850.9, p-value: < 2.22e-16

# Log likelihood: 780.5842 for error model
# ML residual variance (sigma squared): 0.0035464, (sigma: 0.059552)
# Number of observations: 559
# Number of parameters estimated: 4
# AIC: -1553.2, (AIC for lm: -1443.4)



sar_cell <- errorsarlm(en_NA~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
# summary(sar_cell) 

# Call:errorsarlm(formula = en_NA ~ hs, listw = nlw_cell, na.action = na.omit,
    # zero.policy = T)

# Residuals:
       # Min         1Q     Median         3Q        Max
# -0.1032103 -0.0228405 -0.0163772 -0.0072714  0.4821862

# Type: error
# Regions with no neighbours included:
 # 5641 5769 5920 6384 6914 7301 8131 8262 8393 8524 8778 9019 9431 9926 10215 10344 10481 10579 10602 10757 10869 10985 10997 11001 11394 11415 11634 11651 11772 11904 12855 13125 13253 13889 13891
# Coefficients: (asymptotic standard errors)
              # Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.1091772  0.0079214 13.7826 < 2.2e-16
# hs          -0.0478799  0.0116878 -4.0966 4.193e-05

# Lambda: 0.76784, LR test value: 263.31, p-value: < 2.22e-16
# Asymptotic standard error: 0.020962
    # z-value: 36.631, p-value: < 2.22e-16
# Wald statistic: 1341.8, p-value: < 2.22e-16

# Log likelihood: 856.3535 for error model
# ML residual variance (sigma squared): 0.0020348, (sigma: 0.045109)
# Number of observations: 559
# Number of parameters estimated: 4
# AIC: -1704.7, (AIC for lm: -1443.4)


# Now compare how much Variation can be explained
summary(lm_IP_NA)$adj.r.squared # The r_squared of the normal regression: 0.01261041

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.1929991
summary(sar_cell,Nagelkerke=T)$NK # 0.3846225 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_IP_NA) #  780.5842                   724.699, = 0.0009024

LR.sarlm(sar_cell,lm_IP_NA) #  856.3535                   724.6992, p-value = 2.61e-12



### plot correlogram adding sar_IP_NA model
##-------------------------------------------

## residual for sar_IP_i

rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_IP_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_IP_i

rval <- getValues(en_crop_NA) # Create new raster
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
