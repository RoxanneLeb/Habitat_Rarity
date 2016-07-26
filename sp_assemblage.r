
# Sp assemblage #
#################

# 14/07/2016: Roxanne Leberger

# To assemble the dataset of the birds provided by Graeme:
# esh_sp_1.tif = resident birds
# esh_sp_2.tif = breeding area for migrating birds
# esh_sp_3.tif = non-breeding area for migrating birds

# create list

cd ../..

cd net/netmonde1/vol/vol17/leberro-processing/esh_Graeme/

# start R

R

# common_list_1_2_3

library(raster)

list_1 <- Sys.glob(file.path("*_1.tif"))
list_2 <- Sys.glob(file.path("*_2.tif"))
list_3 <- Sys.glob(file.path("*_3.tif"))

strsplit(list_1,"_")  # resident birds
strsplit(list_2,"_")  # breeding area for migrating birds
strsplit(list_3,"_")  # non-breeding area for migrating birds

s_1 <- strsplit(list_1,"_")
list_1_name <- sapply(s_1, function(x) { length(x) <- 3; x[2] })
list_1_name <- as.factor(list_1_name)

s_2 <- strsplit(list_2,"_") #
list_2_name <- sapply(s_2, function(x) { length(x) <- 3; x[2] })
list_2_name <- as.factor(list_2_name)

s_3 <- strsplit(list_3,"_")
list_3_name <- sapply(s_3, function(x) { length(x) <- 3; x[2] })
list_3_name <- as.factor(list_3_name)

#  common_list_1_2, birds with 1 and 2
common_list_1_2 <- intersect(list_1_name,list_2_name)

#  common_list_1_3, birds with 1 and 3
common_list_1_3 <- intersect(list_1_name,list_3_name)

#  common_list_2_3, birds with 2 and 3
common_list_2_3 <- intersect(list_2_name,list_3_name)

# common_list_1_2_3, birds with 1, 2 and 3
common_list_1_2_3 <- Reduce(intersect, list(list_1_name,list_2_name,list_3_name))

length(common_list_1_2) # 961 resident birds with breeding area
length(common_list_1_3) # 1020 resident birds with non breeding area
length(common_list_1_2_3)  # 668 resident birds with breeding and non breeding area
length(common_list_2_3) # 1227

length(list_1) # 9679 resident birds
length(list_2) # 1530 breeding area
length(list_3) # 1582 non breeding area

# sum rasters list_1 + list_2 => list_12, and then sum list_12 to list_3
##~~~~~~~

common_list_1_2[i] 
common_12_list_1 <- paste("esh_", common_list_1_2, "_1.tif", sep = "")
common_12_list_2 <- paste("esh_", common_list_1_2, "_2.tif", sep = "")

HS_forest_realm <- raster('../HS_realm_forest_11var.tif')

system('ls *_12.tif | wc -l')

for(i in 1:length(common_list_1_2)){
rast_list1 <- raster(common_12_list_1[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_12_list_1_proj <- projectRaster(rast_list1, HS_forest_realm)
rast_list2 <- raster(common_12_list_2[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_12_list_2_proj <- projectRaster(rast_list2, HS_forest_realm)

common_12_list_1_proj[is.na(common_12_list_1_proj)] <- 0
common_12_list_2_proj[is.na(common_12_list_2_proj)] <- 0

esh_12_rast <- mosaic(common_12_list_1_proj,common_12_list_2_proj, fun='max')

writeRaster(esh_12_rast,  paste("esh_", common_list_1_2[i], "_12.tif", sep = ""), overwrite=T)
}

# sum the missing rasters list_12+list_3
##~~~~~~

length(unique(common_list_1_2)) # 961
length(unique(union(common_list_1_2,common_list_1_2_3))) # 961 rasters == length(common_list_1_2) => no difference between 2 vectors
length(common_list_1_2_3) # 668 rasters => all rasters are part of comon_list_1_2

# to add rasters with _3 to the _1_2
common_123_list_3 <- paste("esh_", common_list_1_2_3, "_3.tif", sep = "")
common_123_list_12 <- paste("esh_", common_list_1_2_3, "_12.tif", sep = "")

HS_forest_realm <- raster('../HS_realm_forest_11var.tif')

for(i in 1:length(common_123_list_3)){
common_123_list_12_proj <- raster(common_123_list_12[i])
common_123_list_3_rast <- raster(common_123_list_3[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_123_list_3_proj <- projectRaster(common_123_list_3_rast, HS_forest_realm)

common_123_list_3_proj[is.na(common_123_list_3_proj)] <- 0

esh_12_rast <- mosaic(common_123_list_12_proj, common_123_list_3_proj, fun='max')

writeRaster(esh_12_rast,  paste("esh_", common_list_1_2_3[i], "_123.tif", sep = ""))
file.remove(paste("esh_", common_list_1_2_3[i], "_12.tif", sep = ""))
}


system('ls *_123.tif | wc -l')

# sum the missing rasters list_2+list_3
##~~~~~~

exclude_same_123 <- intersect(common_list_2_3, common_list_1_2_3)

length(common_list_2_3 [! common_list_2_3 %in% common_list_1_2_3])
list_2_3 <- common_list_2_3 [! common_list_2_3 %in% common_list_1_2_3]


common_23_list_2 <- paste("esh_", list_2_3, "_2.tif", sep = "")
common_23_list_3 <- paste("esh_", list_2_3, "_3.tif", sep = "")

HS_forest_realm <- raster('../HS_realm_forest_11var.tif')

system('ls *_23.tif | wc -l')

for(i in 1:length(list_2_3)){
rast_list2 <- raster(common_23_list_2[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_23_list_2_proj <- projectRaster(rast_list2, HS_forest_realm)
rast_list3 <- raster(common_23_list_3[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_23_list_3_proj <- projectRaster(rast_list3, HS_forest_realm)

common_23_list_2_proj[is.na(common_23_list_2_proj)] <- 0
common_23_list_3_proj[is.na(common_23_list_3_proj)] <- 0

esh_23_rast <- mosaic(common_23_list_2_proj,common_23_list_3_proj, fun='max')

writeRaster(esh_23_rast,  paste("esh_", list_2_3[i], "_23.tif", sep = ""), overwrite=T)
}

# sum the missing rasters list_1+list_3
##~~~~~~

exclude_same_123 <- intersect(common_list_1_3, common_list_1_2_3)

length(common_list_1_3 [! common_list_1_3 %in% common_list_1_2_3])
list_1_3 <- common_list_1_3 [! common_list_1_3 %in% common_list_1_2_3]

system('ls *_13.tif | wc -l') # 352

common_13_list_1 <- paste("esh_", list_1_3, "_1.tif", sep = "")
common_13_list_3 <- paste("esh_", list_1_3, "_3.tif", sep = "")

HS_forest_realm <- raster('../HS_realm_forest_11var.tif')

for(i in 1:length(list_1_3)){
rast_list1 <- raster(common_13_list_1[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_13_list_1_proj <- projectRaster(rast_list1, HS_forest_realm)
rast_list3 <- raster(common_13_list_3[i], crs=' +proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs')
common_13_list_3_proj <- projectRaster(rast_list3, HS_forest_realm)

common_13_list_1_proj[is.na(common_13_list_1_proj)] <- 0
common_13_list_3_proj[is.na(common_13_list_3_proj)] <- 0

esh_13_rast <- mosaic(common_13_list_1_proj,common_13_list_3_proj, fun='max')

writeRaster(esh_13_rast,  paste("esh_", list_1_3[i], "_13.tif", sep = ""), overwrite=T)
}

# to put all files together without repetition
system('ls *.tif > birds_list_mess.txt')
bird_list <- read.table('birds_list_mess.txt')
dim(bird_list)

list_1 <- Sys.glob(file.path("*_1.tif"))
list_2 <- Sys.glob(file.path("*_2.tif"))
list_3 <- Sys.glob(file.path("*_3.tif"))

list_12 <- Sys.glob(file.path("*_12.tif"))
list_13 <- Sys.glob(file.path("*_13.tif"))
list_23 <- Sys.glob(file.path("*_23.tif"))
list_123 <- Sys.glob(file.path("*_123.tif"))

length(list_1) # 9679 # after cleaning: 8637
length(list_2) # 1530 # after cleaning: 344
length(list_3) # 1582 # after cleaning: 343
length(list_12) # 293 
length(list_13) # 352
length(list_23) # 559
length(list_123) # 668

9679+1530+1582+293+352+559+668 # 14675, 12792 original files

strsplit(list_1,"_")  

s_1 <- strsplit(list_1,"_")
list_1_name <- sapply(s_1, function(x) { length(x) <- 3; x[2] })
list_1_name <- as.factor(list_1_name)

s_2 <- strsplit(list_2,"_")
list_2_name <- sapply(s_2, function(x) { length(x) <- 3; x[2] })
list_2_name <- as.factor(list_2_name)

s_3 <- strsplit(list_3,"_")
list_3_name <- sapply(s_3, function(x) { length(x) <- 3; x[2] })
list_3_name <- as.factor(list_3_name)

s_12 <- strsplit(list_12,"_")
list_12_name <- sapply(s_12, function(x) { length(x) <- 3; x[2] })
list_12_name <- as.factor(list_12_name)

s_13 <- strsplit(list_13,"_")
list_13_name <- sapply(s_13, function(x) { length(x) <- 3; x[2] })
list_13_name <- as.factor(list_13_name)

s_23 <- strsplit(list_23,"_")
list_23_name <- sapply(s_23, function(x) { length(x) <- 3; x[2] })
list_23_name <- as.factor(list_23_name)

s_123 <- strsplit(list_123,"_")
list_123_name <- sapply(s_123, function(x) { length(x) <- 3; x[2] })
list_123_name <- as.factor(list_123_name)

list_name <- c(as.character(list_1_name),as.character(list_2_name),as.character(list_3_name),as.character(list_12_name),as.character(list_13_name),as.character(list_23_name),as.character(list_123_name))
list_name_complete <- c(list_1,list_2,list_3,list_12,list_13,list_23,list_123)

tab_name <- cbind(list_name,list_name_complete)
head(tab_name)

length(unique(list_name)) # 9679
unique_name <- unique(list_name)

tab_name_duplicate <- tab_name[duplicated(tab_name[,1]), ]

head(tab_name_duplicate,20) # 4412

tab_name_duplicate_sort <- tab_name_duplicate[order(tab_name_duplicate[,1]),] 

head(tab_name_duplicate_sort,40)

# removing duplicate files = 4412
# for(i in 1:dim(tab_name_duplicate_sort)[1]){
# file.remove(paste("esh_", tab_name_duplicate_sort[i,][1], "_1.tif", sep = ""))
# file.remove(paste("esh_", tab_name_duplicate_sort[i,][1], "_2.tif", sep = ""))
# file.remove(paste("esh_", tab_name_duplicate_sort[i,][1], "_3.tif", sep = ""))
# }

# to put all files together without repetition
system('ls *.tif > birds_list_order1.txt')
bird_list2 <- read.table('birds_list_order1.txt')
dim(bird_list2) # 10251

list_1 <- Sys.glob(file.path("*_1.tif"))
list_2 <- Sys.glob(file.path("*_2.tif"))
list_3 <- Sys.glob(file.path("*_3.tif"))

list_12 <- Sys.glob(file.path("*_12.tif"))
list_13 <- Sys.glob(file.path("*_13.tif"))
list_23 <- Sys.glob(file.path("*_23.tif"))
list_123 <- Sys.glob(file.path("*_123.tif"))

length(list_1) # 9679 => 8366
length(list_2) # 1530 => 10
length(list_3) # 1582 => 3
length(list_12) # 293 => 293
length(list_13) # 352 => 352
length(list_23) # 559 => 559
length(list_123) # 668 => 668


s_1 <- strsplit(list_1,"_")
list_1_name <- sapply(s_1, function(x) { length(x) <- 3; x[2] })
list_1_name <- as.factor(list_1_name)

s_2 <- strsplit(list_2,"_")
list_2_name <- sapply(s_2, function(x) { length(x) <- 3; x[2] })
list_2_name <- as.factor(list_2_name)

s_3 <- strsplit(list_3,"_")
list_3_name <- sapply(s_3, function(x) { length(x) <- 3; x[2] })
list_3_name <- as.factor(list_3_name)

s_12 <- strsplit(list_12,"_")
list_12_name <- sapply(s_12, function(x) { length(x) <- 3; x[2] })
list_12_name <- as.factor(list_12_name)

s_13 <- strsplit(list_13,"_")
list_13_name <- sapply(s_13, function(x) { length(x) <- 3; x[2] })
list_13_name <- as.factor(list_13_name)

s_23 <- strsplit(list_23,"_")
list_23_name <- sapply(s_23, function(x) { length(x) <- 3; x[2] })
list_23_name <- as.factor(list_23_name)

s_123 <- strsplit(list_123,"_")
list_123_name <- sapply(s_123, function(x) { length(x) <- 3; x[2] })
list_123_name <- as.factor(list_123_name)

list_name <- c(as.character(list_1_name),as.character(list_2_name),as.character(list_3_name),as.character(list_12_name),as.character(list_13_name),as.character(list_23_name),as.character(list_123_name))
list_name_complete <- c(list_1,list_2,list_3,list_12,list_13,list_23,list_123)

tab_name <- cbind(list_name,list_name_complete)
head(tab_name)

length(unique(list_name)) # 10251

unique_name_sort <- tab_name[order(tab_name[,1]),] 

tail(unique_name_sort,150)
