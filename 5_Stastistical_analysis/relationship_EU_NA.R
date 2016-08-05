
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

length(en_EU_mask_0[en_EU_mask_0==0]) # 2240 cells = 0
length(en_EU_mask_0[en_EU_mask_0>0]) # 81 cells are not 0 values ~ 

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



## Autocorrelation model : MORAN I NOT SIGNIFICANT FOR EUROPE !
###--------------------------------


## Correlogram to set the correlation distance for the SAR model
###-------------------------------------------------------------------

library(ncf) # For the Correlogram


# residual raster: lm_EU_NA
##----------------------------------

## residual for lm_EU_NA
rval <- getValues(en_crop_NA) # Create new raster
rval[as.numeric(names(lm_EU_NA$residuals))]<- lm_EU_NA$residuals # replace all data-cells with res value
resid_lm_NA <- HS_EU_crop
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
#Moran I statistic standard deviate = 0.79718, p-value = 0.2127 
#alternative hypothesis: greater
#sample estimates:
#Moran I statistic       Expectation          Variance
#       0.05527065       -0.01612903        0.00802202



# correlogram to set d2 value for the dnearneigh to calculate SARerr
##

system.time(co_lm_NA <- correlog(x,y,z_lm_NA,increment = 500000,resamp = 0, latlon = F,na.rm=T)) # this can take a while.

# plot linear model correlogram
plot(co_lm_NA, xlim=c(100000,5000000))
#plot(0,type="n",col="black",ylab="Moran's I",xlab="lag distance",xlim=c(100000,7000000),ylim=c(-1,1))
abline(h=0,lty="dotted")
lines(co_lm_NA$correlation~co_lm_NA$mean.of.class,col="red",lwd=2)
points(x=co_lm_NA$x.intercept,y=0,pch=19,col="red")
