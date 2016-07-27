
#############################################
## rasterize forested areas for each realm ##
#############################################

# 15/06/16 : Roxanne Leberger

# open R 

cd net/netmonde1/vol/vol17/leberro-processing/HS_realm/HS_cor_var/

library(raster)

pas_1 <- raster('pas/pa_1.tif')
pa1_proj <- projectRaster(pas_1, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa1_proj_bin <- pa1_proj
pa1_proj_bin[pa1_proj_bin>0] <- 1
pa1_proj_bin[pa1_proj_bin!=1] <- 0

pa1_crop <- crop(pa1_proj_bin,extent(pas_1))

writeRaster(pa1_proj_bin, 'pa_1_bin.tif')


pas_2 <- raster('pas/pa_2.tif')
pa2_proj <- projectRaster(pas_2, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa2_proj_bin <- pa2_proj
pa2_proj_bin[pa2_proj_bin>0] <- 1
pa2_proj_bin[pa2_proj_bin!=1] <- 0

pa2_crop <- crop(pa2_proj_bin,extent(pas_2))

writeRaster(pa2_crop, 'pa_2_bin.tif')


pas_4 <- raster('pa_old/pa_4.tif')
pa4_proj <- projectRaster(pas_4, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa4_proj_bin <- pa4_proj
pa4_proj_bin[pa4_proj_bin>0] <- 1
pa4_proj_bin[pa4_proj_bin!=1] <- 0

writeRaster(pa4_proj_bin, 'pa_4_bin.tif')


pas_5 <- raster('pa_old/pa_5.tif')
pa5_proj <- projectRaster(pas_5, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa5_proj_bin <- pa5_proj
pa5_proj_bin[pa5_proj_bin>0] <- 1
pa5_proj_bin[pa5_proj_bin!=1] <- 0

writeRaster(pa5_proj_bin, 'pa_5_bin.tif')


pas_6 <- raster('pa_6.tif')
pa6_proj <- projectRaster(pas_6, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa6_proj_bin <- pa6_proj
pa6_proj_bin[pa6_proj_bin>0] <- 1
pa6_proj_bin[pa6_proj_bin!=1] <- 0

pa2_crop <- crop(pa6_proj_bin,extent(pas_6))

writeRaster(pa6_proj_bin, 'pa_6_bin.tif')


pas_4 <- raster('/pa_4.tif')
pa4_proj <- projectRaster(pas_4, res=5000, crs='+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')

pa4_proj_bin <- pa4_proj
pa4_proj_bin[pa4_proj_bin>0] <- 1
pa4_proj_bin[pa4_proj_bin!=1] <- 0

writeRaster(pa4_proj_bin, 'pa_4_bin.tif')