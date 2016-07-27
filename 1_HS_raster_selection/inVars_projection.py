
################
## Set inVars ##
################

# Set projection, resolution, NA values

## 13/06/16: Roxanne Leberger

# grass70

# ndvimin average
r.in.gdal in=ndvimin.tif out=ndvimin_1km
r.resamp.stats in=ndvimin_1km out=ndvimin_5km method=average
g.region -d
r.mapcalc 'ndvimin_5km_extent = ndvimin_5km'

# write
r.out.gdal in=ndvimin_5km out=ndvimin_5km_NAval.tif nodata=65535

# ndvimax average
##

r.in.gdal in=ndvimax.tif out=ndvimax_1km
r.resamp.stats in=ndvimax_1km out=ndvimax_5km method=average
g.region -d
r.mapcalc 'ndvimax_5km_extent = ndvimax_5km'

# verification
r.info ndvimax_5km_extent
d.mon wx0
d.rast ndvimax_5km_extent

# rename
g.remove rast=ndvimax_5km
g.copy ndvimax_5km_extent,ndvimax_5km
g.remove rast=ndvimax_5km_extent

# write
r.out.gdal in=ndvimax_5km out=ndvimax_5km_NAval.tif nodata=65535


# arid average
##

r.in.gdal in=aridity_moll.tif out=arid_1km
r.resamp.stats in=arid_1km out=arid_5km method=average
g.region -d
r.mapcalc 'arid_5km_extent = arid_5km'

# verification
r.info arid_5km_extent
d.mon wx0
d.rast arid_5km_extent

# rename
g.remove type=rast name=arid_5km -f
g.copy arid_5km_extent,arid_5km --o
g.remove rast=arid_5km_extent -f

# write
r.out.gdal in=arid_5km out=arid_5km_NAval.tif nodata=65535


# ndwi
##

r.in.gdal in=ndwi.tif out=ndwi_1km
r.resamp.stats in=ndwi_1km out=ndwi_5km method=average
g.region -d
r.mapcalc 'ndwi_5km_extent = ndwi_5km'

r.info ndwi_5km_extent
d.mon wx0
d.rast ndwi_5km_extent

# rename
g.remove rast=ndwi_5km
g.copy ndwi_5km_extent,ndwi_5km
g.remove rast=ndwi_5km_extent

# write
r.out.gdal in=ndwi_5km out=ndwi_5km_NAval.tif nodata=255


# pre
r.in.gdal in=pre.tif out=pre_1km
r.resamp.stats in=pre_1km out=pre_5km method=average
g.region -d
r.mapcalc 'pre_5km_extent = pre_5km'

r.info pre_5km_extent
d.mon wx0
d.rast pre_5km_extent

# rename
g.remove type=raster name=pre_5km@leberro -f
g.copy pre_5km_extent,pre_5km
g.remove type=raster name=pre_5km_extent -f

# write
r.out.gdal in=pre_5km out=pre_5km_NAval.tif nodata=65535

# tseas
##

r.in.gdal in=tseason_moll.tif out=tseas_1km
r.resamp.stats in=tseas_1km out=tseas_5km method=average
g.region -d
r.mapcalc 'tseas_5km_extent = tseas_5km'

# rename
g.remove type=raster name=tseas_5km@leberro -f
g.copy tseas_5km_extent,tseas_5km
g.remove type=raster name=tseas_5km_extent -f

r.info tseas_5km
d.mon wx0
d.rast tseas_5km

# write
r.out.gdal in=tseas_5km out=tseas_5km_NAval.tif nodata=65535

# treeP
##

r.in.gdal in=tree.tif out=treeP_1km
r.resamp.stats in=treeP_1km out=treeP_5km method=average
g.region -d
r.mapcalc 'treeP_5km_extent = treeP_5km'

# rename
g.remove type=raster name=treeP_5km@leberro -f
g.copy treeP_5km_extent,treeP_5km
g.remove type=raster name=treeP_5km_extent -f

r.info treeP_5km
d.mon wx0
d.rast treeP_5km

# write
r.out.gdal in=treeP_5km out=treeP_5km_NAval.tif nodata=255

# treeH
r.in.gdal in=tree_height_moll.tif out=treeH_1km
r.resamp.stats in=treeH_1km out=treeH_5km method=average
g.region -d
r.mapcalc 'treeH_5km_extent = treeH_5km'

# rename
g.remove type=raster name=treeH_5km@leberro -f
g.copy treeH_5km_extent,treeH_5km
g.remove type=raster name=treeH_5km_extent -f

r.info treeH_5km
d.mon wx0
d.rast treeH_5km

# write
r.out.gdal in=treeH_5km out=treeH_5km_NAval.tif nodata=65535

# soil
r.in.gdal in=soil_ph_moll.tif out=soil_1km
r.resamp.stats in=soil_1km out=soil_5km method=average
g.region -d
r.mapcalc 'soil_5km_extent = soil_5km'

# rename
g.remove type=raster name=soil_5km@leberro -f
g.copy soil_5km_extent,soil_5km
g.remove type=raster name=soil_5km_extent -f

r.info soil_5km
d.mon wx0
d.rast soil_5km

# write
r.out.gdal in=soil_5km out=soil_5km_NAval.tif nodata=65535

# slope
r.in.gdal in=slope.tif out=slope_1km
r.resamp.stats in=slope_1km out=slope_5km method=average
g.region -d
r.mapcalc 'slope_5km_extent = slope_5km'

# rename
g.remove type=raster name=slope_5km@leberro -f
g.copy slope_5km_extent,slope_5km
g.remove type=raster name=slope_5km_extent -f

r.info slope_5km
d.mon wx0
d.rast slope_5km

# write
r.out.gdal in=slope_5km out=slope_5km_NAval.tif nodata=65535


# temp_max_warm_month_moll
##

r.in.gdal in=temp_max_warm_month_moll.tif out=tmax_warm_M_1km
r.resamp.stats in=tmax_warm_M_1km out=tmax_warm_M_5km method=average
g.region -d
r.mapcalc 'tmax_warm_M_5km_extent = tmax_warm_M_5km'

# verification
r.info tmax_warm_M_5km_extent
d.mon wx0
d.rast tmax_warm_M_5km_extent

# rename
g.remove type=rast name=tmax_warm_M_5km -f
g.copy tmax_warm_M_5km_extent,tmax_warm_M_5km --o
g.remove type=rast name=tmax_warm_M_5km_extent -f

# write
r.out.gdal in=tmax_warm_M_5km out=tmax_warm_M_5km_NAval.tif nodata=65535


# pre_driest_month_moll
r.in.gdal in=pre_driest_month_moll.tif out=preD_M_1km
r.resamp.stats in=preD_M_1km out=preD_M_5km method=average
g.region -d
r.mapcalc 'preD_M_5km_extent = preD_M_5km'

# rename
g.remove type=raster name=preD_M_5km@leberro -f
g.copy preD_M_5km_extent,preD_M_5km
g.remove type=raster name=preD_M_5km_extent -f

r.info preD_M_5km
d.mon wx0
d.rast preD_M_5km

# write
r.out.gdal in=preD_M_5km out=preD_M_5km_NAval.tif nodata=65535


# tree density
##

r.in.gdal in=tree_density_moll.tif out=treeD_1km
r.resamp.stats in=treeD_1km out=treeD_5km method=average
g.region -d
r.mapcalc 'treeD_5km_extent = treeD_5km'

# rename
g.remove type=raster name=treeD_5km@leberro -f
g.copy treeD_5km_extent,treeD_5km
g.remove type=raster name=treeD_5km_extent -f

r.info treeD_5km
d.mon wx0
d.rast treeD_5km

# write
r.out.gdal in=treeD_5km out=treeD_5km_NAval.tif nodata=1.8e+308


