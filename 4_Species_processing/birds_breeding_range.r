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

write.table(birds_1_list_diff,'birds_1_list_toProject.txt',row.names = F, col.names=F)

q()

R
birds <- Sys.glob(file.path("esh_original_files/*_2.tif"))
# write.table(birds,'birds_2_list.txt', col.names=F)
birds_2 <- read.table('birds_2_list.txt')
# write.table(birds_2[,2],'birds_2_list_toProject.txt',row.names = F, col.names=F)

# to copy the selected rasters in another folder

# e:
# cd E:\leberro\My Documents\PhD_Paper_1_globalForest\Database\BirdLife\list_birds

mkdir out_1_2
cd esh_original_files/
# tr -d '\15\32' < birds_1_list_toProject.txt > birds_1_list_toProject2.txt # to remove all carriage returns and Ctrl-z (^Z) characters from a Windows file
cp birds_1_list_toProject.txt birds_1_list_toProject2.txt
cat birds_1_list_toProject2.txt | xargs -I {} cp {} ../out_1_2/
# type birds_1_list_toProject.txt | forfiles /out_1_2_2 copy {} out_1_2_2/ # on windows

cat birds_2_list_toProject.txt | xargs -I {} cp {} ../out_2_2/

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

R

# load forest file
HS_forest_realm <- raster('../HS_realm_forest_11var.tif')


# to sum the 1 in the forested area
#~~~~

ls *_bin.tif > birds_100km_bin.txt 
birds=$(cat "birds_100km_bin.txt")        #the output of 'cat $file' is assigned to the $name variable
echo $birds   

# create mask for forested area
r.mask raster=HS_realm_forest_11var_100km # forest mask

# count cells of each raster present in forest
for file in  $birds  ; do 
r.stats in=$file -c >> range_forest_birds_100km.txt ; # -c = count cell sp range covered by mask, add in file
done ;

r.mask -r

## output :
# 0 : global forests
# 1 : forested bird area => what I want
# * : NA = non forest

# do the same in R

R

library(raster)
library(doParallel)
library(foreach)

HS_realm_forest_100km <- raster('../../Paper1_globalForest/HS_realm/HS_11var/HS_realm_forest_11var_100km.tif')

forest_100km <- HS_realm_forest_100km
forest_100km[!is.na(forest_100km)] <- 1

birds <- list.files(pattern="*_bin.tif$")
n <- length(birds)

range_tab <- data.frame(rep(NA,length(birds)),rep(NA,length(birds)))

# with foreach parallel
##

registerDoParallel(50) 
range_tab_v1 <- foreach(i=1:n, .packages = "raster",.combine="c") %dopar% {
r1 <- raster(birds[i])
r1_forest <- r1+forest_100km
range_tab[,1][i] <- freq(r1_forest,value=2) # select r1_forest=2
}

range_tab_v2 <- foreach(i=1:n, .packages = "raster",.combine="c") %dopar% {
r1 <- raster(birds[i])
range_tab[,2][i] <- freq(r1,value=1) #  select r1_forest=1
}

range_tab2 <- cbind(as.character(birds),as.numeric(range_tab_v1),as.numeric(range_tab_v2))
colnames(range_tab2) <- c('birds','forest_range','all_range')
head(range_tab2,20)

# write.table(range_tab2,'range_tab.txt',row.names=F)

# or loop for
# for(i in 1:n){ # 6704
# r1 <- raster(birds[i])
# r1_forest <- r1+forest_100km
# range_tab[,1][i] <- freq(r1_forest,value=2) # select r1_forest=2
# range_tab[,2][i] <- freq(r1,value=1) #  select r1_forest=1
# }

range_tab <- read.table('range_tab.txt',h=T)
head(range_tab)
range_tab2 <- range_tab[range_tab[,3]>0,] # divided by 0 (non forest sp) generate NA, remove 0
dim(range_tab2) # 9820 sp with range in forest
forest_range <- as.numeric(range_tab2[,1])/as.numeric(range_tab2[,2])
range_tab3 <- cbind(range_tab2,forest_range)
head(range_tab3,30)

hist(range_tab3[,4],100)
quantile(range_tab3[,4])

range_forest_0.3 <-range_tab3[range_tab3[,4]>=0.33,]
dim(range_forest_0.3) # 9125 sp
range_forest_0.5 <-range_tab3[range_tab3[,4]>=0.5,]
dim(range_forest_0.5) # 8651

range_forest_0.75 <-range_tab3[range_tab3[,4]>=0.75,]
dim(range_forest_0.75) # 7278
# write.table(range_forest_0.75, 'range_forest_75%.txt')

range_forest_0.70 <-range_tab3[range_tab3[,4]>=0.70,]
dim(range_forest_0.70) # 7639


# forested range of species
###

R

range_forest <- read.table('range_forest_75%.txt')

head(range_forest,20)

quantile(range_forest[,3])
#  0%  25%  50%  75% 100%
#   1   20   82  289 8384

# selection of species with a range of 20 cells or less
 
small_range_forest <- range_forest[range_forest[,3]<=20,]
dim(small_range_forest) # 2428
head(small_range_forest, 20)

# write.table(small_range_forest,'small_range_forest.txt', row.names=F)

dim(small_range_forest) # 2428
dim(range_forest) # 9679

# sum rasters of range_forest_75% and small_range_forest

small_range_sp_list <- small_range_forest[,1]
forest_sp_list <- range_forest[,1]
head(forest_sp)

# write.table(small_range_sp_list,'small_range_sp_list.txt', col.names=F, row.names=F)
# write.table(forest_sp_list,'forest_sp_list.txt', col.names=F, row.names=F)
