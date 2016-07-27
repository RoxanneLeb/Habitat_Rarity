
grass70

r.in.gdal in=results/555_1.tif out=555_1 -oe
r.in.gdal in=results/555_2.tif out=555_2 -oe
r.in.gdal in=results/555_3.tif out=555_3 -oe
r.in.gdal in=results/555_4.tif out=555_4 -oe
r.in.gdal in=results/555_5.tif out=555_5 -oe
r.in.gdal in=results/555_6.tif out=555_6 -oe


# set all the maps at the same region

g.region -d
g.region -p

r.mapcalc '555_1_reg = 555_1'
r.mapcalc '555_2_reg = 555_2'
r.mapcalc '555_3_reg = 555_3'
r.mapcalc '555_4_reg = 555_4'
r.mapcalc '555_5_reg = 555_5'
r.mapcalc '555_6_reg = 555_6'

d.mon wx0
d.rast 555_1_reg

# replace NA values by 0 (to make sum)

r.null map=555_1_reg null=0
r.null map=555_2_reg null=0
r.null map=555_3_reg null=0
r.null map=555_4_reg null=0
r.null map=555_5_reg null=0
r.null map=555_6_reg null=0

# add maps to create one

r.out.gdal in=555_1_reg out=555_1_reg.tif
r.out.gdal in=555_2_reg out=555_2_reg.tif
r.out.gdal in=555_3_reg out=555_3_reg.tif
r.out.gdal in=555_4_reg out=555_4_reg.tif
r.out.gdal in=555_5_reg out=555_5_reg.tif
r.out.gdal in=555_6_reg out=555_6_reg.tif