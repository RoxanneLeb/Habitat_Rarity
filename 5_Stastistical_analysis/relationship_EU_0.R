
## EU Relationship
###--------------------------
####-------------------------


## load raster
###---------------------------

HS_realm2<- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var/HS_realm_11var_redifined_100km.tif')
en <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/sp_maps/endemic_birds_perc.tif')

# define study area: EUstraly
EU <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/EU_01_realm_redifined2.tif') # to define EU realm
plot(EU)

EU_100km <- projectRaster(EU, HS_realm2)
EU_100km[EU_100km>=0]<-0

HS_EU_mask <- mask(HS_realm2, EU_100km)
en_EU_mask <- mask(en, EU_100km)

# put 0 in forest without small range birds
en_EU_mask_0 <- en_EU_mask
en_EU_mask_0[is.na(en_EU_mask) & !is.na(HS_EU_mask)] <- 0

# crop for faster computation
##

e <- c(-5000000,180000000,1000000,10000000)

HS_EU_crop <- crop(HS_EU_mask,e)
plot(HS_EU_crop)

en_crop_0 <- crop(en_EU_mask_0,e)
plot(en_crop_0)

en_crop_NA <- crop(en_EU_mask,e)
plot(en_crop_NA)


## Linear model with 0 and NA
###------------------------------

table(HS_EU_mask[HS_EU_mask==0]) # ? cells = 0
length(en_EU_mask[en_EU_mask>0]) # 81 cells are not 0 values ~ 

lm_EU_NA <- lm(values(en_crop_NA)~values(HS_EU_crop),na.action=na.omit)
summary(lm_EU_NA)
# Multiple R-squared:  0.04236,  Adjusted R-squared:  0.03024 
# F-statistic: 3.494 on 1 and 79 DF,  p-value: 0.06528 NS
plot(values(en_crop_NA)~values(HS_EU_crop))

lm_EU_0 <- lm(values(en_crop_0)~values(HS_EU_crop),na.action=na.omit)
summary(lm_EU_0)
# Multiple R-squared:  0.01558,  Adjusted R-squared:  0.01515 
# F-statistic:  36.7 on 1 and 2319 DF,  p-value: 1.607e-09 ***

plot(lm_EU_0)

plot(values(en_crop_0)~values(HS_EU_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the Paleartic realm')

abline(0.022389,-0.028916 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0019287,-0.0024501,col='orange') # when 0 kept = regression in all forested area
points(x=HS_EU_mask[en_EU_mask_0==0],y=en_EU_mask_0[en_EU_mask_0==0],col='orange')



## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_EU_0
##

## residual for lm_EU_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_EU_0$residuals))]<- lm_EU_0$residuals # replace all data-cells with res value
resid_lm_0 <- HS_EU_crop
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
# Moran I statistic standard deviate = 9.5493, p-value < 2.2e-16
# alternative hypothesis: greater
# sample estimates:
# Moran I statistic       Expectation          Variance
#      0.1228819489     -0.0004372540      0.0001667715
 


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

co_lm_0$x.intercept # 2111541
co_lm_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_EU_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
en_0 = getValues(en_crop_0)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=2111541,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_EU_crop)


# spatial error SAR

sar_i <- errorsarlm(en_0~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_i, na.action = na.omit,
    # zero.policy = T)

# Residuals:
        # Min          1Q      Median          3Q         Max
# -0.00517169 -0.00061362 -0.00034275 -0.00020153  0.21833231

# Type: error
# Coefficients: (asymptotic standard errors)
               # Estimate  Std. Error z value Pr(>|z|)
# (Intercept)  0.01132919  0.00584151  1.9394  0.05245
# hs          -0.00057760  0.00048254 -1.1970  0.23131

# Lambda: 0.97853, LR test value: 107.88, p-value: < 2.22e-16
# Asymptotic standard error: 0.014281
    # z-value: 68.518, p-value: < 2.22e-16
# Wald statistic: 4694.7, p-value: < 2.22e-16

# Log likelihood: 8564.001 for error model
# ML residual variance (sigma squared): 3.6333e-05, (sigma: 0.0060277)
# Number of observations: 2321
# Number of parameters estimated: 4
# AIC: -17120, (AIC for lm: -17014)


sar_cell <- errorsarlm(en_0~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_cell, na.action = na.omit,
    # zero.policy = T)

# Residuals:
        # Min          1Q      Median          3Q         Max
# -0.00615076 -0.00111183 -0.00031062  0.00015644  0.22025013

# Type: error
# Regions with no neighbours included:
 # 2470 2928 3986 4668 4671 4787 4900 4942 5019 5241 5251 5363 5422 5592 6088 6170 6408 6542 6547 6634 7501 7736 8196 8222 8345 8361 8602 8704 9066 9089 9517 10103 10892
# Coefficients: (asymptotic standard errors)
               # Estimate  Std. Error z value  Pr(>|z|)
# (Intercept)  0.00231728  0.00033040  7.0137 2.322e-12
# hs          -0.00287574  0.00051689 -5.5635 2.644e-08

# Lambda: 0.29025, LR test value: 81.423, p-value: < 2.22e-16
# Asymptotic standard error: 0.02487
    # z-value: 11.671, p-value: < 2.22e-16
# Wald statistic: 136.21, p-value: < 2.22e-16

# Log likelihood: 8550.773 for error model
# ML residual variance (sigma squared): 3.6023e-05, (sigma: 0.0060019)
# Number of observations: 2321
# Number of parameters estimated: 4
# AIC: -17094, (AIC for lm: -17014)



# Now compare how much Variation can be explained
summary(lm_EU_0)$adj.r.squared # The r_squared of the normal regression: 0.01515271

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.06028549 => best
summary(sar_cell,Nagelkerke=T)$NK # 0.04951284


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_EU_0) #8564.001                  8510.062, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_EU_0) # 8550.773                   8510.062, p-value < 2.2e-16



### plot correlogram adding sar_EU_0 model

## residual for sar_EU_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_EU_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_EU_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_cell$residuals))]<- sar_cell$residuals # replace all data-cells with res value
resid_sar_cell <- HS_EU_crop
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
