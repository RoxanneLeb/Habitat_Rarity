##########################
# forest birds selection #
##########################

# 26/07/2016: Roxanne Leberger

# To select species with 33% of their breeding area covered by forest on a 100km resolution basis

# work on gosling (server)

cd ../..

cd net/netmonde1/vol/vol17/leberro-processing/esh_Graeme/


###--- optional section -------------

# To select _1.tif and _2.tif sp wich are not processed yet to change projection
###

# starting R

R

library(raster)


# select _1.tif birds which were not reproject yet

birds <- Sys.glob(file.path("esh_original_files/*_1.tif"))
length(birds) # 9679 of _1.tif

# write.table(birds,'birds_1_list.txt')
birds_1 <- read.table('birds_1_list.txt')
birds <- as.data.frame(birds_1[,2])

# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/BirdLife/list_birds/')

birds_1_2 <- Sys.glob(file.path("out_1_2/*_1.tif"))
# write.table(birds_1_2,'birds_1_2_list.txt')
birds_1_2_proj <- read.table('birds_1_2_list.txt')


# difference between 2 lists

head(birds)
head(birds_1_2_proj)

birds_1_num <- strsplit(as.character(birds[,1]),"_")
birds_1_cat <- sapply(birds_1_num, function(x) {length(x) <- 3; x[2] })
birds_1_num <- as.factor(birds_1_cat)

birds_1_num_proj <- strsplit(as.character(birds_1_2_proj[,1]),"_")
birds_1_cat_proj <- sapply(birds_1_num_proj, function(x) {length(x) <- 6; x[5] })
birds_1_num_proj <- as.factor(birds_1_cat_proj)

diff1 <- setdiff(birds_1_num, birds_1_num_proj)
length(diff1) #
length(birds_1_num)-length(birds_1_num_proj)

birds_1_list_diff <- paste("esh_", diff1, "_1.tif", sep = "")
birds_1_list_diff

write.table(birds_1_list_diff,'birds_1_list_toProject.txt',row.names = F)

q()

# to copy the selected rasters in another folder

# e:
# cd E:\leberro\My Documents\PhD_Paper_1_globalForest\Database\BirdLife\list_birds

mkdir out_1_2_2
cd esh_original_files/
tr -d '\15\32' < birds_1_list_toProject.txt > birds_1_list_toProject2.txt # to remove all carriage returns and Ctrl-z (^Z) characters from a Windows file
cp birds_1_list_toProject.txt birds_1_list_toProject2.txt
cat birds_1_list_toProject2.txt | xargs -I {} cp {} ../out_1_2_2/
# type birds_1_list_toProject.txt | forfiles /out_1_2_2 copy {} out_1_2_2/ # on windows

###--- optional section END ------------------


# ArcGIS for projection: use Iterator
###

# Change resolution with grass70
###

grass70

# loop: resolution change : 100km

cd out_1_2/
birds=$(cat "out_1_diff.txt")        #the output of 'cat $file' is assigned to the $name variable
echo $birds   

g.region res=100000 # set resolution at 100km

for file in  $birds; do 
	r.in.gdal in=$file'.tif' out=$file --o
	r.resamp.stats input=$file  output=$file'_100km_diff'  method=average --o
	r.null map=$file'_100km_diff' null=0 --o
done;

# select species with 33% of their breeding area covered by forest
####

# range of these species

range <- rep(NA,length(birds))

for(i in 1:length(birds)){
r1 <- raster(birds[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
r1_cpt <- r1
r1_cpt[r1_cpt > 0] <- 1
range[i] <- sum(values(!is.na(r1_cpt))) # 258846
}

write.table(range,'range_breeding_birds.txt')

# load forest file
HS_forest_realm <- raster('../HS_realm_forest_11var.tif')

# remove areas of 0 cells
range_tab2 <- range_tab[range_tab[,2] != 0,]
dim(range_tab2)