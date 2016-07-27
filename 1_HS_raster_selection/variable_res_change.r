# R- transform variable raster #
################################

## 18/05/2016 - Roxanne Leberger

# using R on server gosling 
# connect on hanks 22 with putty ssh -X leberro@0.0.0.0 -p 27 password: ugis

library(raster)

setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/var_5km/')

# agregate 1km -> 5km resolution

pre <- raster('pre.tif')

# writeRaster(pre_5km, 'pre_5km_moll.tif')

ndvimin <- raster('ndvimin.tif')
ndvimin_5km <- aggregate(ndvimin, 5, fun=mean)
ndvimin_5km_2 <- projectRaster(ndvimin_5km, pre_5km)
# writeRaster(ndvimin_5km_2, 'ndvimin_5km_moll.tif')

ndvimax <- raster('ndvimax.tif')
ndvimax_5km <- aggregate(ndvimax, 5, fun=mean)
ndvimax_5km_2 <- projectRaster(ndvimax_5km, pre_5km)
# writeRaster(ndvimax_5km_2, 'ndvimax_5km_moll.tif', overwrite=T)

slope <- raster('slope.tif')
slope_5km <- aggregate(slope, 5, fun=mean)
# writeRaster(slope_5km, 'slope_5km_moll.tif')

ndwi <- raster('ndwi.tif')
ndwi_5km <- aggregate(ndwi, 5, fun=mean)
# writeRaster(ndwi_5km, 'ndwi_5km_moll.tif')

tree <- raster('tree.tif')
tree_5km <- aggregate(tree, 5, fun=mean)
# writeRaster(tree_5km, 'tree_5km_moll.tif')

herb <- raster('herb.tif')
herb_5km <- aggregate(herb, 5, fun=mean)
# writeRaster(herb_5km, 'herb_5km_moll.tif')

# to put all rasters at the same resolution
##~~~~~~~~

pre_5km <- raster('pre_5km_moll.tif') # as reference for projection

arid <- raster('aridity_moll.tif')
arid_5km_1 <- aggregate(arid, 5, fun=mean)
arid_5km_2 <- projectRaster(arid_5km_1, pre_5km)
arid_5km_2
arid_mask <- mask(arid_5km_2,pre_5km)
plot(arid_mask)
# writeRaster(arid_mask, 'arid5km_moll.tif', overwrite=T)

soilPH <- raster('soil_ph_moll.tif')
soilPH_5km_1 <- aggregate(soilPH, 5, fun=mean)
soilPH_5km_2 <- projectRaster(soilPH, pre_5km)
soilPH_5km_2
soilPH_mask <- mask(soilPH_5km_2,pre_5km)
plot(soilPH_mask)
# writeRaster(soilPH_mask, 'soilPH5km_moll.tif', overwrite=T)

tmean <- raster('tmean_moll.tif')
tmean_5km_1 <- aggregate(tmean, 5, fun=mean)
tmean_5km_2 <- projectRaster(tmean_5km_1, pre_5km)
tmean_5km_2
tempM_mask <- mask(tmean_5km_2,pre_5km)
# writeRaster(tempM_mask, 'tmean5km_moll.tif', overwrite=T)

tseason <- raster('tseason_moll.tif')
tseason_5km_1 <- aggregate(tseason, 5, fun=mean)
tseason_5km_2 <- projectRaster(tseason_5km_1, pre_5km)
tseason_5km_2
tseason_mask <- mask(tseason_5km_2,pre_5km)
plot(tseason_mask)
# writeRaster(tseason_mask, 'tseason5km_moll.tif', overwrite=T)

tsums0 <- raster('tsums0_moll.tif')
tsums0_5km_1 <- aggregate(tsums0, 5, fun=mean)
tsums0_5km_2 <- projectRaster(tsums0_5km_1, pre_5km)
tsums0_5km_2
tsums0_mask <- mask(tsums0_5km_2,pre_5km)
plot(tsums0_mask)
# writeRaster(tsums0_mask, 'tsums05km_moll.tif', overwrite=T)

tree_height <- raster('tree_height_moll.tif')
tree_height_5km_1 <- aggregate(tree_height, 5, fun=mean)
tree_height_5km_2 <- projectRaster(tree_height_5km_1, pre_5km)
tree_height_5km_2
treeH_mask <- mask(tree_height_5km_2,pre_5km)
plot(treeH_mask)
# writeRaster(treeH_mask, 'treeH5km_moll.tif', overwrite=T)

petseason <- raster('petseason_moll.tif')
petseason_5km_1 <- aggregate(petseason, 5, fun=mean)
petseason_5km_2 <- projectRaster(petseason_5km_1, pre_5km)
petseason_5km_2
petseason_mask <- mask(petseason_5km_2,pre_5km)
plot(petseason_mask)
# writeRaster(petseason_mask, 'petseason5km_moll.tif', overwrite=T)


treed <- raster('Tree_density_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
treed_5km_1 <- aggregate(treed, 5, fun=mean)
treed_5km_2 <- projectRaster(treed_5km_1, pre_5km)
treed_5km_2
treed_mask <- mask(treed_5km_2,pre_5km)
plot(treed_mask)
# writeRaster(treed_mask, 'treeD_5km_moll.tif', overwrite=T)


setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Variables/To_be_transformed')

preD_month <- raster('pre_driest_month_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
preD_month_5km_1 <- aggregate(preD_month, 5, fun=mean)
preD_month_5km_2 <- projectRaster(preD_month_5km_1, pre_5km)
preD_month_5km_2
preD_month_mask <- mask(preD_month_5km_2,pre_5km)
plot(preD_month_mask)
# writeRaster(preD_month_mask, 'preD_month_5km_moll.tif', overwrite=T)


preD_quarter <- raster('pre_driest_quarter_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
preD_quarter_5km_1 <- aggregate(preD_quarter, 5, fun=mean)
preD_quarter_5km_2 <- projectRaster(preD_quarter_5km_1, pre_5km)
preD_quarter_5km_2
preD_quarter_mask <- mask(preD_quarter_5km_2,pre_5km)
plot(preD_quarter_mask)
# writeRaster(preD_quarter_mask, 'preD_quarter_5km_moll.tif', overwrite=T)


twarm_quarter <- raster('temp_mean_warmest_quarter_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
twarm_quarter_5km_1 <- aggregate(twarm_quarter, 5, fun=mean)
twarm_quarter_5km_2 <- projectRaster(twarm_quarter_5km_1, pre_5km)
twarm_quarter_5km_2
twarm_quarter_mask <- mask(twarm_quarter_5km_2,pre_5km)
plot(twarm_quarter_mask)
# writeRaster(twarm_quarter_mask, 'twarm_quarter_5km_moll.tif', overwrite=T)


trange <- raster('temp_annual_range_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
trange_5km_1 <- aggregate(trange, 5, fun=mean)
trange_5km_2 <- projectRaster(trange_5km_1, pre_5km)
trange_5km_2
trange_mask <- mask(trange_5km_2,pre_5km)
plot(trange_mask)
# writeRaster(trange_mask, 'trange_5km_moll.tif', overwrite=T)


tmin_cold_month <- raster('temp_min_cold_month_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
tmin_cold_month_5km_1 <- aggregate(tmin_cold_month, 5, fun=mean)
tmin_cold_month_5km_2 <- projectRaster(tmin_cold_month_5km_1, pre_5km)
tmin_cold_month_5km_2
tmin_cold_month_mask <- mask(tmin_cold_month_5km_2,pre_5km)
plot(tmin_cold_month_mask)
# writeRaster(tmin_cold_month_mask, 'tmin_cold_month_5km_moll.tif', overwrite=T)


tmax_warm_month <- raster('temp_max_warm_month_moll.tif')
pre_5km <- raster('pre_5km_NAval.tif')
tmax_warm_month_5km_1 <- aggregate(tmax_warm_month, 5, fun=mean)
tmax_warm_month_5km_2 <- projectRaster(tmax_warm_month_5km_1, pre_5km)
tmax_warm_month_5km_2
tmax_warm_month_mask <- mask(tmax_warm_month_5km_2,pre_5km)
plot(tmax_warm_month_mask)
# writeRaster(tmax_warm_month_mask, 'tmax_warm_month_5km_moll.tif', overwrite=T)


