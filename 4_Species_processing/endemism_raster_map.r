
################################
# generate endemism raster map #
################################

# 29/07/2016 : Roxanne Leberger


# work on gosling (server)

cd ../..

cd net/netmonde1/vol/vol17/leberro-processing/esh_Graeme/

grass70

# sh rmapcalc_sum_forest_birds_list.txt
r.info forest_birds
d.mon wx1
d.rast forest_birds
r.null map=forest_birds setnull=0 # replace 0 to NA

# to export the map well defined
r.mask HS_realm_forest_11var_100km 
r.out.gdal in=forest_birds out=forest_birds.tif --o

# sh rmapcalc_small_range_sp_sum.txt
r.info small_range_birds
d.mon wx0
d.rast small_range_birds
r.null map=small_range_birds setnull=0
r.out.gdal in=small_range_birds out=forest_birds_small_range.tif

r.mask -r

# open R
R

small_range <- raster('forest_birds_small_range.tif')
birds <- raster('forest_birds.tif')

endemic_perc <- small_range/birds

plot(endemic_perc)

# writeRaster(endemic_perc,'endemic_birds_perc.tif')

# small range / sp
r.mapcalc 'endemic_birds_perc = small_range_birds/forest_birds'
r.info endemic_birds_perc
d.mon wx2
d.rast endemic_birds_perc

##--- END  endemic_birds_perc map ---

### endemism richness map

R

range_forest <- read.table('range_forest_75%.txt')
