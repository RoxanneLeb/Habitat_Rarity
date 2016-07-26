######################################################
# Changing raster resolution and produce binary maps #
######################################################

# 26/07/2016 : Roxanne Leberger

# using gosling (server)

cd ../..

cd net/netmonde1/vol/vol17/leberro-processing/esh_Graeme/

# build list to work on in grass
cd out_1/ 
ls *_1_out.tif > out_1_list.txt
cd ..
cd out_12/ 
ls *_12_out.tif > out_12_list.txt
cd ..
cd out_13/ 
ls *_13_out.tif > out_13_list.txt
cd ..
cd out_23/ 
ls *_23_out.tif > out_23_list.txt
cd ..
cd out_123/ 
ls *_123_out.tif > out_123_list.txt
cd ..


# import list changing resolution

grass70

cd out_1/
birds=$(cat "out_1_list.txt")        #the output of 'cat $file' is assigned to the $name variable
echo $birds   

birds2=$(cat "out_1_list2.txt")    
echo $birds2

# changing resolution loop => doesn't work
# g.region res=100000

# for file in  $birds2; do 
	# r.in.gdal in=$file'_out.tif' out=$file --o
	# r.resamp.stats input=$file  output=rast_100km  method=average --o
	# r.null map=rast_100km null=0 --o
	# r.mapcalc 'rast_100km_2 = rast_100km*1000' --o
	# r.reclass input=rast_100km_2  output=$file'_100km_01' rules=r_reclass_rules_01.txt --o
# done; # not good results... split the loop?

# loop: resolution
g.region res=100000

for file in  $birds2; do 
	r.in.gdal in=$file'_out.tif' out=$file --o
	r.resamp.stats input=$file  output=$file'_100km'  method=average --o
	r.null map=$file'_100km' null=0 --o
done;

d.mon wx0
d.rast mol_esh_22678135_1_100km # ok
d.rast mol_esh_22678135_1 # ok 

# generate binary rasters (1/0 => presence/abs)
##

# r.mapcalc and r.reclass used by sh namefile.txt

paste out_1_list.txt out_1_list2.txt > rmapcalc_out_1_list.txt
head rmapcalc_out_1_list.txt
gedit rmapcalc_out_1_list.txt # transform the names to have the command lines

# sh rmapcalc_out_1_list.txt # execute the file

r.info mol_esh_22678135_1_100km2

# sh rreclass_out_1_list.txt  # execute the file

d.mon wx0 
d.rast mol_esh_22678135_1_100km_bin



