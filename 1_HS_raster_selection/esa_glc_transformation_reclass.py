
# projectRaster GRASS
r.in.gdal in=esa_glc_moll_5km_resample.tif out=esa_glc_moll_5km_resample -oe
g.region -d
g.region -p
r.mapcalc 'esa_glc_5km4 = esa_glc_moll_5km_resample'
r.info esa_glc_5km4
r.out.gdal in=esa_glc_5km4 out=esa_glc_moll5km.tif

# select forest class GRASS
r.in.gdal in=esa_glc_moll5km.tif out=esa_glc_moll5km -oe --o

r.info esa_glc_moll5km

d.mon wx0

r.reclass in=esa_glc_moll5km out=esa_forest_moll5km rules=esa_forest_reclass.txt --o
r.info esa_forest_moll5km

r.out.gdal in=esa_forest_moll5km out=esa_forest_moll5km.tif --o