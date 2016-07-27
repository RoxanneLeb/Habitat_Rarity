library(raster)

getwd()

setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_all_var')

# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-arid')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimax')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimin')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-pre')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-preD_M')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-slope')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-soilPH')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tmax_warm_M')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeD')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeH')
# setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tseas')


reg_555_1 <- raster('555_1_reg.tif')
reg_555_2 <- raster('555_2_reg.tif')
reg_555_3 <- raster('555_3_reg.tif')
reg_555_4 <- raster('555_4_reg.tif')
reg_555_5 <- raster('555_5_reg.tif')
reg_555_6 <- raster('555_6_reg.tif')


# mosaic - need to have the same origin -
reg_555 <- mosaic(reg_555_1,reg_555_2,reg_555_3,reg_555_4,reg_555_5,reg_555_6, fun=max)

# replace 0 by NA
reg_555_NA <- reg_555
reg_555_NA[reg_555 == 0] <- NA

# writeRaster(reg_555_NA, 'HS_10var_tmax_warm_forest.tif', overwrite=T)

# with category

## realm 1

# b to NA
reg_555_1_NA <- reg_555_1
reg_555_1_NA[reg_555_1 == 0] <- NA

# cat
reg_555_1_NA_cat <- reg_555_1_NA
reg_555_1_NA_cat[reg_555_1_NA_cat < quantile(reg_555_1_NA,na.rm = T)[2]] <- 25
reg_555_1_NA_cat[reg_555_1_NA_cat < quantile(reg_555_1_NA,na.rm = T)[3]] <- 50
reg_555_1_NA_cat[reg_555_1_NA_cat < quantile(reg_555_1_NA,na.rm = T)[4]] <- 75
reg_555_1_NA_cat[reg_555_1_NA_cat < 25] <- 100

plot(reg_555_1_NA_cat)


## realm 2

# b to NA
reg_555_2_NA <- reg_555_2
reg_555_2_NA[reg_555_2 == 0] <- NA

# cat
reg_555_2_NA_cat <- reg_555_2_NA
reg_555_2_NA_cat[reg_555_2_NA_cat < quantile(reg_555_2_NA,na.rm = T)[2]] <- 25
reg_555_2_NA_cat[reg_555_2_NA_cat < quantile(reg_555_2_NA,na.rm = T)[3]] <- 50
reg_555_2_NA_cat[reg_555_2_NA_cat < quantile(reg_555_2_NA,na.rm = T)[4]] <- 75
reg_555_2_NA_cat[reg_555_2_NA_cat < 25] <- 100

plot(reg_555_2_NA_cat)


## realm 3

# b to NA
reg_555_3_NA <- reg_555_3
reg_555_3_NA[reg_555_3 == 0] <- NA

# cat
reg_555_3_NA_cat <- reg_555_3_NA
reg_555_3_NA_cat[reg_555_3_NA_cat < quantile(reg_555_3_NA,na.rm = T)[2]] <- 25
reg_555_3_NA_cat[reg_555_3_NA_cat < quantile(reg_555_3_NA,na.rm = T)[3]] <- 50
reg_555_3_NA_cat[reg_555_3_NA_cat < quantile(reg_555_3_NA,na.rm = T)[4]] <- 75
reg_555_3_NA_cat[reg_555_3_NA_cat < 25] <- 100

plot(reg_555_3_NA_cat)


## realm 4

# b to NA
reg_555_4_NA <- reg_555_4
reg_555_4_NA[reg_555_4 == 0] <- NA

# cat
reg_555_4_NA_cat <- reg_555_4_NA
reg_555_4_NA_cat[reg_555_4_NA_cat < quantile(reg_555_4_NA,na.rm = T)[2]] <- 25
reg_555_4_NA_cat[reg_555_4_NA_cat < quantile(reg_555_4_NA,na.rm = T)[3]] <- 50
reg_555_4_NA_cat[reg_555_4_NA_cat < quantile(reg_555_4_NA,na.rm = T)[4]] <- 75
reg_555_4_NA_cat[reg_555_4_NA_cat < 25] <- 100

plot(reg_555_4_NA_cat)



## realm 5

# b to NA
reg_555_5_NA <- reg_555_5
reg_555_5_NA[reg_555_5 == 0] <- NA

# cat
reg_555_5_NA_cat <- reg_555_5_NA
reg_555_5_NA_cat[reg_555_5_NA_cat < quantile(reg_555_5_NA,na.rm = T)[2]] <- 25
reg_555_5_NA_cat[reg_555_5_NA_cat < quantile(reg_555_5_NA,na.rm = T)[3]] <- 50
reg_555_5_NA_cat[reg_555_5_NA_cat < quantile(reg_555_5_NA,na.rm = T)[4]] <- 75
reg_555_5_NA_cat[reg_555_5_NA_cat < 25] <- 100

plot(reg_555_5_NA_cat)


## realm 6

# b to NA
reg_555_6_NA <- reg_555_6
reg_555_6_NA[reg_555_6 == 0] <- NA

# cat
reg_555_6_NA_cat <- reg_555_6_NA
reg_555_6_NA_cat[reg_555_6_NA_cat < quantile(reg_555_6_NA,na.rm = T)[2]] <- 25
reg_555_6_NA_cat[reg_555_6_NA_cat < quantile(reg_555_6_NA,na.rm = T)[3]] <- 50
reg_555_6_NA_cat[reg_555_6_NA_cat < quantile(reg_555_6_NA,na.rm = T)[4]] <- 75
reg_555_6_NA_cat[reg_555_6_NA_cat < 25] <- 100

table(values(reg_555_6_NA_cat))
plot(reg_555_6_NA_cat)



# mosaic - need to have the same origin -
reg_555 <- mosaic(reg_555_1_NA_cat,reg_555_2_NA_cat,reg_555_3_NA_cat,reg_555_4_NA_cat,reg_555_5_NA_cat,reg_555_6_NA_cat, fun=max)


# just Africa
Af <- raster('555_3.tif')

# b to NA
Af_NA <- Af
Af_NA[Af == 0] <- NA

# writeRaster(Af_NA, 'HS_Af_11var.tif')

# cat
Af_NA_cat <- Af_NA
Af_NA_cat[Af_NA_cat < quantile(Af_NA,na.rm = T)[2]] <- 25
Af_NA_cat[Af_NA_cat < quantile(Af_NA,na.rm = T)[3]] <- 50
Af_NA_cat[Af_NA_cat < quantile(Af_NA,na.rm = T)[4]] <- 75
Af_NA_cat[Af_NA_cat < 25] <- 100

table(values(Af_NA_cat))
plot(Af_NA_cat)

# writeRaster(Af_NA_cat, 'HS_Af_11var_cat.tif')


## redefine realms 4 and 1

reg_555_12 <- raster('555_122_reg.tif')
reg_555_42 <- raster('555_422_reg.tif')

# replace 0 by NA
reg_555_12_NA <- reg_555_12
reg_555_12_NA[reg_555_12 == 0] <- NA

# replace 0 by NA
reg_555_42_NA <- reg_555_42
reg_555_42_NA[reg_555_42 == 0] <- NA
plot(reg_555_42_NA)

# writeRaster(reg_555_12_NA, 'HS_PA_redifine_11var.tif', overwrite=T)
# writeRaster(reg_555_42_NA, 'HS_IP_redifine_11var2.tif', overwrite=T)


# mosaic - need to have the same origin -
reg_555 <- mosaic(reg_555_12,reg_555_2,reg_555_3,reg_555_42,reg_555_5,reg_555_6, fun=max)

# replace 0 by NA
reg_555_NA <- reg_555
reg_555_NA[reg_555 == 0] <- NA

HS_realm_forest_11var<-raster('HS_realm_forest_11var.tif')
reg_555_NA2<- projectRaster(reg_555_NA,HS_realm_forest_11var)

# writeRaster(reg_555_NA, 'HS_11var_realm_redifined_forest2.tif', overwrite=T)
