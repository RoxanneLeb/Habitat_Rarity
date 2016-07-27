# R- correlation input variable  #
##################################

## 18/05/2016 - Roxanne Leberger

# using R on server gosling 
# connect on hanks 22 with putty ssh -X leberro@0.0.0.0 -p 27 password: ugis

library(raster)

setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/')

pre <- raster('pre_5km_moll.tif') 
ndvimin <- raster('ndvimin_5km_moll.tif')
ndvimax <- raster('ndvimax_5km_moll.tif')
slope <- raster('slope_5km_moll.tif')
ndwi <- raster('ndwi_5km_moll.tif')
#treeP <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/tree_5km_moll.tif')
#grassP <- raster('herb_5km_moll.tif')
arid <- raster('arid5km_moll.tif')
soil <- raster('soilPH5km_moll.tif')
tmean <- raster('tmean_5km_moll.tif')
tseason <- raster('tseason_5km_moll.tif')
tsums <- raster('tsums0_5km_moll.tif')
treeH <- raster('treeH_5km_moll.tif')
#treeH2 <- projectRaster(treeH, tsums)
# writeRaster(treeH, 'treeH5km_moll.tif')
petseason <- raster('petseason_5km_moll.tif')
treeD <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/treeD_5km_moll.tif')
preseason <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/pre_seas_5km_moll.tif')

preD_Mt <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/preD_month_5km_moll.tif')
preD_Qt <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/preD_quarter_5km_moll.tif')
twarm_Qt <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/twarm_quarter_5km_moll.tif')
tmin_cold_Mt <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/tmin_cold_month_5km_moll.tif')
tmax_warm_Mt <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/tmax_warm_month_5km_moll.tif')
trange <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/trange_5km_moll.tif')

## not servor
# pre <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/pre_5km_moll.tif')
# ndvimin <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/ndvimin_5km_moll.tif')
# ndvimax <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/ndvimax_5km_moll.tif')
# slope <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/slope_5km_moll.tif')
# ndwi <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/ndwi_5km_moll.tif')
# treeP <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/tree_5km_moll.tif')
# grassP <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/herb_5km_moll.tif')
# arid <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/arid_5km_moll.tif')
# soil <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/soilPH_5km_moll.tif')
# tmean <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/tmean_5km_moll.tif')
# tseason <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/tseason_5km_moll.tif')
# tsums <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/tsums0_5km_moll.tif')
# treeH <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/treeH_5km_moll.tif')
# petseason <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/petseason_5km_moll.tif')

## crop the forest part of each raster
forest_01 <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/ESA_GLC/forest_cover_01.tif') # = eco_555

forest <- forest_01
forest[values(forest_01) == 0] <- NA

# for each variable


tmeanf <- tmean
tmeanf[is.na(values(forest))] <- NA
plot(tmeanf)

tseasonf <- tseason
tseasonf[is.na(values(forest))] <- NA
plot(tseasonf)

tsumsf <- tsums
tsumsf[is.na(values(forest))] <- NA
plot(tsumsf)

pref <- pre
pref[is.na(values(forest))] <- NA
plot(pref)

preDf <- preD
preDf[is.na(values(forest))] <- NA
plot(preDf)

preseasonf <- preseason
preseasonf[is.na(values(forest))] <- NA
plot(preseasonf)

petseasonf <- petseason
petseasonf[is.na(values(forest))] <- NA
plot(petseasonf)

ndviminf <- ndvimin
ndviminf[is.na(values(forest))] <- NA
plot(ndviminf)

ndvimaxf <- ndvimax
ndvimaxf[is.na(values(forest))] <- NA
plot(ndvimaxf)

ndwif <- ndwi
ndwif[is.na(values(forest))]<- NA
plot(ndwif)

treeDf <- treeD
treeDf[is.na(values(forest))]<- NA
plot(treeDf)

treeHf <- treeH
treeHf[is.na(values(forest))]<- NA
plot(treeHf)

treePf <- treeP
treePf[is.na(values(forest))] <- NA
plot(treePf)

#grassPf <- grassP
#grassPf[is.na(values(forest))] <- NA
#plot(grassPf)

aridf<- arid
aridf[is.na(values(forest))] <- NA
plot(aridf)

soilf <- soil
soilf[is.na(values(forest))] <- NA
plot(soilf)

slopef <- slope
slopef[is.na(values(forest))] <- NA
plot(slopef)

preD_Mtf <- preD_Mt
preD_Mtf[is.na(values(forest))] <- NA
plot(preD_Mtf)

preD_Qtf <- preD_Qt
preD_Qtf[is.na(values(forest))] <- NA
plot(preD_Qtf)

twarm_Qtf <- twarm_Qt 
twarm_Qtf[is.na(values(forest))] <- NA
plot(twarm_Qtf)

tmin_cold_Mtf <- tmin_cold_Mt 
tmin_cold_Mtf[is.na(values(forest))] <- NA
plot(tmin_cold_Mtf)

tmax_warm_Mtf <- tmax_warm_Mt 
tmax_warm_Mtf[is.na(values(forest))] <- NA
plot(tmax_warm_Mtf)

trangef <- trange
trangef[is.na(values(forest))] <- NA
plot(trangef)


# writeRaster(tmeanf, 'tmean_forest.tif')
# writeRaster(tseasonf, 'tseason_forest.tif')
# writeRaster(tsumsf, 'tsums0_forest.tif')
# writeRaster(pref, 'pre_forest.tif')
# writeRaster(preseasonf, 'preseas_forest.tif')
# writeRaster(ndviminf, 'ndvimin_forest.tif')
# writeRaster(ndvimaxf, 'ndvimax_forest.tif')
# writeRaster(ndwif, 'ndwi_forest.tif')
# writeRaster(treeHf, 'treeH_forest.tif')
# writeRaster(treeDf, 'treeD_forest.tif')
# writeRaster(aridf, 'arid_forest.tif')
# writeRaster(soilf, 'soil_forest.tif')
# writeRaster(slopef, 'slope_forest.tif')
# writeRaster(preD_Mtf, 'preD_Mt_forest.tif')
# writeRaster(preD_Qtf, 'preD_Qt_forest.tif')
# writeRaster(twarm_Qtf, 'twarm_Qt_forest.tif')
# writeRaster(tmin_cold_Mtf, 'tmin_cold_Mt_forest.tif')
# writeRaster(tmax_warm_Mtf, 'tmax_warm_Mt_forest.tif')
# writeRaster(trangef, 'trange_forest.tif')

tmeanf<- raster('tmean_forest.tif')
tseasonf<- raster('tseason_forest.tif')
tsumsf<- raster('tsums0_forest.tif')
tmax_warm_Mtf<- raster('tmax_warm_Mt_forest.tif')
pref<- raster('pre_forest.tif')
preD_Mtf<- raster('preD_Mt_forest.tif')
preseasonf<- raster('preseas_forest.tif')
ndviminf<- raster('ndvimin_forest.tif')
ndvimaxf<- raster('ndvimax_forest.tif')
#ndwif<- raster('ndwi_forest.tif')
treeHf<- raster('treeH_forest.tif')
treeDf<- raster('treeD_forest.tif')
aridf<- raster('arid_forest.tif')
soilf<- raster('soil_forest.tif')
slopef<- raster('slope_forest.tif')

# create dataframe

tmean.df <- data.frame(values(tmeanf))
tseason.df <- data.frame(values(tseasonf))
tsums.df <- data.frame(values(tsumsf))
pre.df <- data.frame(values(pref))
petseason.df <- data.frame(values(petseasonf))
ndvimin.df <- data.frame(values(ndviminf))
ndvimax.df <- data.frame(values(ndvimaxf))
ndwi.df <- data.frame(values(ndwif))
treeH.df <- data.frame(values(treeHf))
treeD.df <- data.frame(values(treeDf))
treeP.df <- data.frame(values(treePf))
arid.df <- data.frame(values(aridf))
soil.df <- data.frame(values(soilf))
slope.df <- data.frame(values(slopef))
preD.df <- data.frame(values(preDf))
preseason.df <- data.frame(values(preseasonf))

preD_Mt.df <- data.frame(values(preD_Mtf))
preD_Qt.df <- data.frame(values(preD_Qtf))
twarm_Qt.df <- data.frame(values(twarm_Qtf))
tmin_cold_Mt.df <- data.frame(values(tmin_cold_Mtf))
tmax_warm_Mt.df <- data.frame(values(tmax_warm_Mtf))
trange.df <- data.frame(values(trangef))

# all_var

var.df <- cbind(tmean.df,
tseason.df,
tsums.df,
pre.df,
petseason.df,
ndvimin.df,
ndvimax.df,
ndwi.df,
treeH.df,
treeD.df,
#grassP.df,
arid.df,
soil.df,
slope.df,
preD.df,
preseason.df)

head(var.df)
summary(var.df)

# selected_var

var.df <- cbind(#tmean.df,
                tseason.df,
                #tsums.df,
                pre.df,
                preD_Mt.df,
                arid.df,
                #preD_Qt.df,
                #twarm_Qt.df,
                #tmin_cold_Mt.df,
                tmax_warm_Mt.df,
                #trange.df,
                #preseason.df,
                #petseason.df,
                ndvimin.df,
                ndvimax.df,
                #ndwi.df,
                treeH.df,
                #treeP.df,
                treeD.df,
                #grassP.df,
                soil.df,
                slope.df)

head(var.df)
summary(var.df)



# correlation test
##

library(corrplot)

var.df2 <- na.omit(var.df)
head(var.df2)
var_cor <- cor(var.df2)

# half plot
var_cor_number <- var_cor
png('var_cor.png')
corrplot(var_cor_number, method="number", type='lower',cl.length=15, number.cex=0.7)
dev.off()

pdf('var_cor.pdf')
corrplot(var_cor_number, method="number", type='lower',cl.length=15, number.cex=0.7)
dev.off()

# mixed plot
var_cor_mixed <- var_cor
png('var_cor.png')
corrplot.mixed(var_cor_mixed, cl.length=15, number.cex=0.7)
dev.off()

pdf('var_cor.pdf')
corrplot.mixed(var_cor_mixed)
dev.off()

#plot(values(bio.df),values(pre.df))
cor.test(ndwi.df[!is.na(ndwi.df)],epr.df[!is.na(ndwi.df)]) #

lm_bio <- lm(values(bio.df)~values(epr.df))
summary(lm_bioepr)


