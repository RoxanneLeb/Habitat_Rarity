
#####################
# Build Birds table #
#####################

# 28/07/2016: Roxanne Leberger

# work on gosling (server)

cd ../..

cd net/netmonde1/vol/vol17/leberro-processing/esh_Graeme/
# select _1.tif birds which were not reproject yet


# create tables with R
#######

R

birds_1_all <- Sys.glob(file.path("esh_original_files/*_1.tif"))
length(birds_1_all) # 9679 of _1.tif
# write.table(birds_1_all,'birds_1_list.txt', col.names=F)

birds_2_all <- Sys.glob(file.path("esh_original_files/*_2.tif"))
length(birds_2_all) # 1530 of _2.tif
# write.table(birds_2_all,'birds_2_list.txt', col.names=F)

birds_1_all <- read.table('birds_1_all_list.txt')
birds_2_all <- read.table('birds_2_all_list.txt')

# build _12.tif list
###

birds_1_num <- strsplit(as.character(birds_1_all[,2]),"_")
birds_1_cat <- sapply(birds_1_num, function(x) {length(x) <- 3; x[2] })
birds_1_num2 <- as.factor(birds_1_cat)

birds_2_num <- strsplit(as.character(birds_2_all[,2]),"_")
birds_2_cat <- sapply(birds_2_num, function(x) {length(x) <- 3; x[2] })
birds_2_num2 <- as.factor(birds_2_cat)

birds_12 <- intersect(birds_1_num2,birds_2_num2)
length(birds_12)

birds_12_list <- paste("esh_", birds_12, "_12.tif", sep = "")
birds_12_list

# write.table(birds_12_list,'birds_12_list.txt',row.names = F, col.names=F)


# build _1_only.tif and _2_only.tif
###

birds_12 <- read.table('birds_12_list.txt')

birds_12_num <- strsplit(as.character(birds_12[,1]),"_")
birds_12_cat <- sapply(birds_12_num, function(x) {length(x) <- 3; x[2] })
birds_12_num2 <- as.factor(birds_12_cat)


# 1_only
#---
birds_1_strict <- setdiff(birds_1_num2,birds_12_num2)
length(birds_1_strict)

birds_1_only_list <- paste("esh_", birds_1_strict, "_1_only.tif", sep = "")
birds_1_only_list

# write.table(birds_1_only_list,'birds_1_only_list.txt',row.names = F, col.names=F)


# 2_only
#---
birds_2_strict <- setdiff(birds_2_num2,birds_12_num2)
length(birds_2_strict)

birds_2_only_list <- paste("esh_", birds_2_strict, "_2_only.tif", sep = "")
birds_2_only_list

# write.table(birds_2_only_list,'birds_2_only_list.txt',row.names = F, col.names=F)

# Create folder moll_1_2_12 with .tif to work with
###

# mkdir moll_1_2_12

cd moll_1
cat birds_1_only_extract_list_out.txt | xargs -I {} cp {} ../moll_1_2_12

# rename _out.tif .tif *_out.tif 

cd moll_2/
cat birds_2_only_extract_list.txt | xargs -I {} cp {} ../moll_1_2_12/

cd moll_1/
cat birds_1_12_extract_list.txt | xargs -I {} cp {} ../moll_12/
cat birds_1_12_extract_list_out.txt | xargs -I {} cp {} ../moll_12/

cd moll_2/
cat birds_2_12_extract_list.txt | xargs -I {} cp {} ../moll_12/

cd moll_12/

# merge _1 and _2 in _12.tif : mosaic
###

R

# for _12 birds

setwd('../moll_12/')
library(raster)

# system('pwd')
# system('ls *.tif > birds_list_12.txt')

# list_1 <- Sys.glob(file.path("*_1.tif"))
# list_2 <- Sys.glob(file.path("*_2.tif"))

s_1 <- strsplit(list_1,"_")
list_1_name <- sapply(s_1, function(x) { length(x) <- 3; x[3] })
list_1_name <- as.factor(list_1_name)

# s_2 <- strsplit(list_2,"_") #
# list_2_name <- sapply(s_2, function(x) { length(x) <- 3; x[3] })
# list_2_name <- as.factor(list_2_name)

# birds_12_list_1 <- paste("mol_esh_", list_1_name, "_1.tif", sep = "")
# birds_12_list_1
# birds_12_list_2 <- paste("mol_esh_", list_2_name, "_2.tif", sep = "")
# birds_12_list_2

HS_forest_realm <- raster('../../Paper1_globalForest/HS_realm/HS_11var/HS_realm_forest_11var.tif')


# parallel computation R
###

#install.packages('foreach')
#install.packages('doParallel')

library(foreach)
library(doParallel)

detectCores() # to have total number of cores (80)
#cl <- makeCluster(5) # choose how many cores
#registerDoParallel(cl)
getDoParWorkers()

#list1 <- list.files(pattern="^.*\\.tif$")
list1 <- list.files(pattern="*_1.tif$")
n <- length(list_1)
list2 <- list.files(pattern="*_2.tif$")

registerDoParallel(30) #or as many cores you want

foreach(i=1:n, .packages = "raster") %dopar% {
birds_12_1_rast <- raster(list_1[i])
birds_12_1_rast_proj <- projectRaster(birds_12_1_rast, HS_forest_realm) # to set a common origin
birds_12_2_rast <- raster(list_2[i])
birds_12_2_rast_proj <- projectRaster(birds_12_2_rast, HS_forest_realm) 
esh_12_rast <- mosaic(birds_12_1_rast_proj,birds_12_2_rast_proj, fun='min')
file_rast <- paste("mol_esh_", list_1_name[i], "_12.tif", sep = "")
writeRaster(esh_12_rast, filename=file_rast)
}

# s_1 <- strsplit(list_1,"_")
# list_1_name <- sapply(s_1, function(x) { length(x) <- 3; x[3] })
# list_1_name <- as.factor(list_1_name)

list12 <- list.files(pattern="*_12.tif$")
s_12 <- strsplit(list12,"_")
list_12_name <- sapply(s_12, function(x) { length(x) <- 3; x[3] })
list_12_name <- as.factor(list_12_name)

setdiff(list_1_name,list_12_name)
# [1] "22698042" "22698049" "22698055"


# normal loop to do the same thing 
for(i in 1:length(list_1)){
birds_12_1_rast <- raster(list_1[i])
birds_12_1_rast_proj <- projectRaster(birds_12_1_rast, HS_forest_realm) # to set a common origin
birds_12_2_rast <- raster(list_2[i])
birds_12_2_rast_proj <- projectRaster(birds_12_2_rast, HS_forest_realm) 
esh_12_rast <- mosaic(birds_12_1_rast_proj,birds_12_2_rast_proj, fun='min')
writeRaster(esh_12_rast,  paste("mol_esh_", list_1_name[i], "_12.tif", sep = ""), overwrite=T)
}


