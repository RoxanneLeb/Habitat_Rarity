
## SAm RelationshSAm
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

length(HS_SAm_mask[HS_SAm_mask==0]) # ? cells = 0
length(en_SAm_mask[en_SAm_mask>0]) # 821 cells are not 0 values ~ 

lm_SAm_NA <- lm(values(en_crop_NA)~values(HS_SAm_crop),na.action=na.omit)
summary(lm_SAm_NA)
# MultSAmle R-squared:  0.07774,  Adjusted R-squared:  0.07661 
# F-statistic: 69.03 on 1 and 819 DF,  p-value: 3.996e-16 ***
plot(values(en_crop_NA)~values(HS_SAm_crop))

lm_SAm_0 <- lm(values(en_crop_0)~values(HS_SAm_crop),na.action=na.omit)
summary(lm_SAm_0)
# MultSAmle R-squared:  0.07034,  Adjusted R-squared:  0.06975 
# F-statistic: 118.9 on 1 and 1572 DF,  p-value: < 2.2e-16

plot(lm_SAm_0)

plot(values(en_crop_0)~values(HS_SAm_crop), xlab='HS values', ylab='percentage of endemism', main='Percentage of endemic birds vs HS values in the SAmrican realm')

abline(0.022389,-0.028916 ,col='red') # when 0 removed = regression just inside sites where small species occur
abline(0.0019287,-0.0024501,col='orange') # when 0 kept = regression in all forested area
points(x=HS_SAm_mask[en_SAm_mask_0==0],y=en_SAm_mask_0[en_SAm_mask_0==0],col='orange')




## Autocorrelation model
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###~~~~~~~~

library(ncf) # For the Correlogram

# Generate an Residual Raster from the linear model
###


# lm_SAm_0
##

## residual for lm_SAm_0
rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(lm_SAm_0$residuals))]<- lm_SAm_0$residuals # replace all data-cells with res value
resid_lm_0 <- HS_SAm_crop
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
#Moran I statistic standard deviate = 46.909, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#Moran I statistic       Expectation          Variance
#     0.9099934948     -0.0006410256      0.0003768559



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

co_lm_0$x.intercept # 3706430
co_lm_0$correlation 


# SAR model
###~~~~~~~~~~~~~~~

# SARerr_SAm_0
#--

## create a list of neighbors

library(spdep)

x = xFromCell(en_crop_0,1:ncell(en_crop_0))
y = yFromCell(en_crop_0,1:ncell(en_crop_0))
en_0 = getValues(en_crop_0)

nb_i <- dnearneigh(cbind(x,y),d1=0,d2=3706430,longlat=F)
nlw_i <- nb2listw(nb_i,style="W",zero.policy=T)

nb_cell <- dnearneigh(cbind(x,y),d1=0,d2=100000,longlat=F)
nlw_cell<- nb2listw(nb_cell,style="W",zero.policy=T)

hs <- values(HS_SAm_crop)


# spatial error SAR

sar_i <- errorsarlm(en_0~hs,listw=nlw_i,na.action=na.omit,zero.policy=T)
summary(sar_i) 

# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_i, na.action = na.omit,
    # zero.policy = T)

# Residuals:
       # Min         1Q     Median         3Q        Max
# -0.0410842 -0.0118754 -0.0027703  0.0015432  0.3151662

# Type: error
# Coefficients: (asymptotic standard errors)
             # Estimate Std. Error  z value Pr(>|z|)
# (Intercept)  0.028516   0.043513   0.6553   0.5123
# hs          -0.031986   0.002224 -14.3822   <2e-16

# Lambda: 0.98416, LR test value: 171.91, p-value: < 2.22e-16
# Asymptotic standard error: 0.011046
    # z-value: 89.095, p-value: < 2.22e-16
# Wald statistic: 7937.9, p-value: < 2.22e-16

# Log likelihood: 3432.212 for error model
# ML residual variance (sigma squared): 0.00074413, (sigma: 0.027279)
# Number of observations: 1574
# Number of parameters estimated: 4
# AIC: -6856.4, (AIC for lm: -6686.5)


sar_cell <- errorsarlm(en_0~hs,listw=nlw_cell,na.action=na.omit,zero.policy=T)
summary(sar_cell) 
# Call:errorsarlm(formula = en_0 ~ hs, listw = nlw_cell, na.action = na.omit,
    # zero.policy = T)

# Residuals:
       # Min         1Q     Median         3Q        Max
# -0.1014335 -0.0021706 -0.0010816  0.0003363  0.1024677

# Type: error
# Regions with no neighbours included:
 # 2466 2643 2793 3228 3393 3529 9025 9627 13839 13991 14137 14289 14446
# Coefficients: (asymptotic standard errors)
              # Estimate Std. Error z value  Pr(>|z|)
# (Intercept)  0.0138810  0.0020605  6.7366 1.621e-11
# hs          -0.0040462  0.0021759 -1.8595   0.06295

# Lambda: 0.9017, LR test value: 2812, p-value: < 2.22e-16
# Asymptotic standard error: 0.0090705
    # z-value: 99.41, p-value: < 2.22e-16
# Wald statistic: 9882.4, p-value: < 2.22e-16

# Log likelihood: 4752.278 for error model
# ML residual variance (sigma squared): 9.8235e-05, (sigma: 0.0099114)
# Number of observations: 1574
# Number of parameters estimated: 4
# AIC: -9496.6, (AIC for lm: -6686.5)



# Now compare how much Variation can be explained
summary(lm_SAm_0)$adj.r.squared # The r_squared of the normal regression:  0.06974728

summary(sar_i,Nagelkerke=T)$NK # Nagelkerkes pseudo r_square of the SAR: 0.1665266
summary(sar_cell,Nagelkerke=T)$NK # 0.8442485 => best


# Finally do a likelihood ratio test
LR.sarlm(sar_i,lm_SAm_0) # 3432.212                   3346.257, p-value < 2.2e-16

LR.sarlm(sar_cell,lm_SAm_0) #   4752.278                   3346.257, p-value < 2.2e-16



### plot correlogram adding sar_SAm_0 model

## residual for sar_SAm_i

rval <- getValues(en_crop_0) # Create new raster
rval[as.numeric(names(sar_i$residuals))]<- sar_i$residuals # replace all data-cells with res value
resid_sar_i <- HS_SAm_crop
values(resid_sar_i) <- rval
rm(rval) #replace our values in this new raster
names(resid_sar_i) <- "Residuals"


## residual for sar_SAm_i

rval <- getValues(en_crop_0) # Create new raster
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
