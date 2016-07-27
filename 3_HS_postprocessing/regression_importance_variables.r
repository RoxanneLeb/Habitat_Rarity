
setwd('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var')

arid <- raster('HS_10var_arid_forest.tif')
ndvimax <- raster('HS_10var_ndvimax_forest.tif')
ndvimin <- raster('HS_10var_ndvimin_forest.tif')
pre <- raster('HS_10var_pre_forest.tif')
preD <- raster('HS_10var_preD_forest.tif')
slope <- raster('HS_10var_slope_forest.tif')
soilPH <- raster('HS_10var_soilPH_forest.tif')
treeH <- raster('HS_10var_treeH_forest.tif')
tmax <- raster('HS_10var_tmax_warm_forest.tif')
treeD <- raster('HS_10var_treeD_forest.tif')
tseas <- raster('HS_10var_tseas_forest.tif')

all <- raster('../HS_all_var/HS_realm_forest_11var.tif')

lm_arid <- lm(values(all) ~ values(arid))
lm_ndvimax <- lm(values(all) ~ values(ndvimax))
lm_ndvimin <- lm(values(all) ~ values(ndvimin))
lm_pre <- lm(values(all) ~ values(pre))
lm_preD <- lm(values(all) ~ values(preD))
lm_slope <- lm(values(all) ~ values(slope))
lm_soilPH <- lm(values(all) ~ values(soilPH))
lm_treeH <- lm(values(all) ~ values(treeH))
lm_tmax <- lm(values(all) ~ values(tmax))
lm_treeD <- lm(values(all) ~ values(treeD))
lm_tseas<- lm(values(all) ~ values(tseas))

#lm_all <- lm(values(all) ~ values(all))

summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)
summary(lm_all)


plot(20:30,20:30)
abline(0.6327746, 1, col='red')
abline(2.320e-02, 9.981e-01, col='orange') # aridity
abline(0.0132467, 0.9865888, col='darkgreen') # ndvimax
abline(0.0210322, 0.9694387, col='lightgreen') # ndvimin
abline(0.0195405, 0.9934185, col='blue') # pre
abline(0.0185480, 0.9714902, col='lightblue') # preD
abline(0.0173888, 0.9623213, col='grey') # slope
abline(0.0452086, 0.9094410, col='brown') # soilpH
abline(0.0146288, 0.9760817, col='turquoise') # treeH
abline(0.0148165, 0.9900391, col='gold') # tmax
abline(0.0314444, 0.9685557, col='purple') # treeD
abline(0.0096706, 0.9858011, col='pink') # tseas


## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               2.320e-02,
               0.0132467,
               0.0210322,
               0.0195405,
               0.0185480,
               0.0173888,
               0.0452086,
               0.0146288,
               0.0148165,
               0.0314444,
               0.0096706)
pente <- c(1,
           9.981e-01, 
           0.9865888, 
           0.9694387, 
           0.9934185, 
           0.9714902, 
           0.9623213, 
           0.9094410, 
           0.9760817,
           0.9900391,
           0.9685557,
           0.9858011)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))


# lm all
lm_all <- lm(values(all) ~ values(arid)+values(ndvimax))
summary(lm_all)

# arid_dif <- sqrt( ( values(all) - values(all)-values(arid) ) ^2) # the same than lm_all <- lm(values(all) ~ values(arid)+values(ndvimax))
# ndvimax_dif <- sqrt( ( values(all) - values(all)-values(ndvimax) ) ^2) 
# lm_all <- lm(values(all) ~ arid_dif+ndvimax_dif)
# summary(lm_all)
# plot(0:10,0:10)
# abline(-2.182e-02, 1, col='red')
# abline(-2.182e-02, 1.309e-01, col='orange') # aridity
# abline(-2.182e-02, 1.111e-01, col='darkgreen') # ndvimax


lm_all <- lm(values(all) ~ values(arid)+values(ndvimax)+values(ndvimin)+values(pre)+values(preD)+values(slope)+values(soilPH)+values(treeH)+values(tmax)+values(treeD)+values(tseas))
summary(lm_all)
# 
plot(0:10,0:10)
abline(-2.182e-02, 1, col='red')
abline(-2.182e-02, 1.309e-01, col='orange') # aridity
abline(-2.182e-02, 1.111e-01, col='darkgreen') # ndvimax
abline(-2.182e-02, 1.005e-01, col='lightgreen') # ndvimin
abline(-2.182e-02, 1.354e-02, col='blue') # pre
abline(-2.182e-02, 1.293e-01, col='lightblue') # preD
abline(-2.182e-02, 1.256e-01, col='grey') # slope
abline(-2.182e-02, 1.520e-02, col='brown') # soilpH
abline(-2.182e-02, 1.022e-01, col='turquoise') # treeH
abline(-2.182e-02, 8.622e-02, col='gold') # tmax
abline(-2.182e-02, 8.691e-02, col='purple') # treeD
abline(-2.182e-02, 1.404e-01, col='pink') # tseas

library(MASS)
stepAIC(lm_all)
# Start:  AIC=-13721644
# values(all) ~ values(arid) + values(ndvimax) + values(ndvimin) + 
#   values(pre) + values(preD) + values(slope) + values(soilPH) + 
#   values(treeH) + values(tmax) + values(treeD) + values(tseas)
# 
# Df Sum of Sq     RSS       AIC
# <none>                          866.91 -13721644
# - values(pre)      1      0.65  867.55 -13720306
# - values(soilPH)   1      3.66  870.56 -13714085
# - values(arid)     1     41.86  908.77 -13636906
# - values(tmax)     1     58.96  925.86 -13603425
# - values(slope)    1    133.51 1000.42 -13464262
# - values(treeD)    1    152.51 1019.41 -13430466
# - values(treeH)    1    162.54 1029.45 -13412864
# - values(ndvimax)  1    163.68 1030.59 -13410883
# - values(ndvimin)  1    200.60 1067.51 -13347636
# - values(tseas)    1    245.94 1112.85 -13272896
# - values(preD)     1    343.64 1210.55 -13121694
# 
# Call:
#   lm(formula = values(all) ~ values(arid) + values(ndvimax) + values(ndvimin) + 
#        values(pre) + values(preD) + values(slope) + values(soilPH) + 
#        values(treeH) + values(tmax) + values(treeD) + values(tseas))

# ANOVA a 11 facteurs



## For each realm
#################

## Realm 01 - Paleartic

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')

arid_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-arid/555_122.tif')
ndvimax_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimax/555_122.tif')
ndvimin_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimin/555_122.tif')
pre_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-pre/555_122.tif')
preD_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-preD_M/555_122.tif')
slope_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-slope/555_122.tif')
soilPH_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-soilPH/555_122.tif')
treeH_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeH/555_122.tif')
tmax_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tmax_warm_M/555_122.tif')
treeD_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeD/555_122.tif')
tseas_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tseas/555_122.tif')

all_PA <- all
all_PA[arid_PA==0] <- NA


lm_arid <- lm(values(all_PA) ~ values(arid_PA))
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))

#lm_all <- lm(values(all_PA) ~ values(all_PA))

summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)

## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.0960566,
               0.0845580,
               0.1096316,
               0.0903861,
               0.0677418,
               0.1025507,
               0.1212633,
               0.0790185,
               0.0836816,
               0.1154345,
               0.0877283)
pente <- c(1,
           0.9339155, 
           0.9363405, 
           0.9066696, 
           0.9408497, 
           0.9633844, 
           0.8914616, 
           0.8569256, 
           0.9345162,
           0.9321225,
           0.9006449,
           0.9302662)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))




## Realm 02 - Australia

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')
PA <- realm
PA[realm !=2] <- NA
plot(PA)

all_PA <- all
all_PA[is.na(PA)] <- NA

arid_PA <- arid
arid_PA[is.na(PA)] <- NA

ndvimax_PA <- ndvimax
ndvimax_PA[is.na(PA)] <- NA

ndvimin_PA <- ndvimin
ndvimin_PA[is.na(PA)] <- NA

pre_PA <- pre
pre_PA[is.na(PA)] <- NA

preD_PA <- preD
preD_PA[is.na(PA)] <- NA

slope_PA <- slope
slope_PA[is.na(PA)] <- NA

soilPH_PA <- soilPH
soilPH_PA[is.na(PA)] <- NA

treeH_PA <- treeH
treeH_PA[is.na(PA)] <- NA

tmax_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tmax_warm_M/555_2_reg.tif')
tmax_PA[is.na(PA)] <- NA

treeD_PA <- treeD
treeD_PA[is.na(PA)] <- NA

tseas_PA <- tseas
tseas_PA[is.na(PA)] <- NA


lm_arid <- lm(values(all_PA) ~ values(arid_PA), na.rm=T)
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))

#lm_all <- lm(values(all_PA) ~ values(all_PA))

summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)
summary(lm_all)


plot(20:30,20:30)
abline(0, 1, col='red')
abline(0.0178054, 0.9875300, col='orange') # aridity
abline(0.021395, 0.956945, col='darkgreen') # ndvimax
abline(0.011277 , 0.995203, col='lightgreen') # ndvimin
abline(0.0201596, 0.9867218, col='blue') # pre
abline(0.017927, 0.982458, col='lightblue') # preD
abline(0.019930, 0.954311, col='grey') # slope
abline(0.052343, 0.886877, col='brown') # soilpH
abline(0.0209366, 0.9717735, col='turquoise') # treeH
abline(0.0229335, 0.9667299, col='gold') # tmax
abline(0.035286, 0.954420, col='purple') # treeD
abline(0.0052653, 0.9862477, col='pink') # tseas


## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.0178054,
               0.021395,
               0.011277,
               0.0201596,
               0.017927,
               0.019930,
               0.052343,
               0.0209366,
               0.0229335,
               0.035286,
               0.0052653)
pente <- c(1,
           0.9875300, 
           0.956945, 
           0.995203, 
           0.9867218, 
           0.982458, 
           0.954311, 
           0.886877, 
           0.9717735,
           0.9667299,
           0.954420,
           0.9862477)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))



## Realm 04 - Africa

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')
PA <- realm
PA[realm != 4] <- NA
plot(PA)

all_PA <- all
all_PA[is.na(PA)] <- NA

arid_PA <- arid
arid_PA[is.na(PA)] <- NA

ndvimax_PA <- ndvimax
ndvimax_PA[is.na(PA)] <- NA

ndvimin_PA <- ndvimin
ndvimin_PA[is.na(PA)] <- NA

pre_PA <- pre
pre_PA[is.na(PA)] <- NA

preD_PA <- preD
preD_PA[is.na(PA)] <- NA

slope_PA <- slope
slope_PA[is.na(PA)] <- NA

soilPH_PA <- soilPH
soilPH_PA[is.na(PA)] <- NA

treeH_PA <- treeH
treeH_PA[is.na(PA)] <- NA

tmax_PA <- tmax
tmax_PA[is.na(PA)] <- NA

treeD_PA <- treeD
treeD_PA[is.na(PA)] <- NA

tseas_PA <- tseas
tseas_PA[is.na(PA)] <- NA



lm_arid <- lm(values(all_PA) ~ values(arid_PA))
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))


summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)


plot(20:30,20:30)
abline(0, 1, col='red')
abline(0.0273005, 0.9844471, col='orange') # aridity
abline(0.0251593, 0.9572420, col='darkgreen') # ndvimax
abline(0.0326679, 0.9430746, col='lightgreen') # ndvimin
abline(0.0268446, 0.9760094, col='blue') # pre
abline(0.0521201, 0.9126687, col='lightblue') # preD
abline(0.0262129, 0.9542288, col='grey') # slope
abline(0.0585646, 0.8763324, col='brown') # soilpH
abline(0.0295239, 0.9506707, col='turquoise') # treeH
abline(0.0301159, 0.9499788, col='gold') # tmax
abline(0.0375960, 0.9821117, col='purple') # treeD
abline(0.0309227, 0.9458264, col='pink') # tseas


## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.0273005,
               0.0251593,
               0.0326679,
               0.0268446,
               0.0521201,
               0.0262129,
               0.0585646,
               0.0295239,
               0.0301159,
               0.0375960,
               0.0309227
               )
pente <- c(1,
           0.9844471, 
           0.9572420, 
           0.9430746, 
           0.9760094, 
           0.9126687, 
           0.9542288, 
           0.8763324, 
           0.9506707,
           0.9499788,
           0.9821117,
           0.9458264
           )

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))



## Realm 06 - South Am

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')
PA <- realm
PA[realm != 6] <- NA
plot(PA)

all_PA <- all
all_PA[is.na(PA)] <- NA

arid_PA <- arid
arid_PA[is.na(PA)] <- NA

ndvimax_PA <- ndvimax
ndvimax_PA[is.na(PA)] <- NA

ndvimin_PA <- ndvimin
ndvimin_PA[is.na(PA)] <- NA

pre_PA <- pre
pre_PA[is.na(PA)] <- NA

preD_PA <- preD
preD_PA[is.na(PA)] <- NA

slope_PA <- slope
slope_PA[is.na(PA)] <- NA

soilPH_PA <- soilPH
soilPH_PA[is.na(PA)] <- NA

treeH_PA <- treeH
treeH_PA[is.na(PA)] <- NA

tmax_PA <- tmax
tmax_PA[is.na(PA)] <- NA

treeD_PA <- treeD
treeD_PA[is.na(PA)] <- NA

tseas_PA <- tseas
tseas_PA[is.na(PA)] <- NA



lm_arid <- lm(values(all_PA) ~ values(arid_PA))
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))


summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)


## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.0382216,
               0.0139740,
               0.0099522,
               0.0292841,
               0.0269034,
               0.0054095,
               0.0376577,
               0.0171329,
               0.0090155,
               0.0242073,
               -0.0116484)
pente <- c(1,
           0.9923436, 
           0.9937953, 
           0.9770683, 
           0.9822599, 
           0.9401097, 
           0.9740402, 
           0.9167314, 
           0.9698806,
           1.0127300,
           0.9939435,
           1.0192954)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11


p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))


# Realm 05 -

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')
PA <- realm
PA[realm != 7] <- NA
plot(PA)

all_PA <- all
all_PA[is.na(PA)] <- NA

arid_PA <- arid
arid_PA[is.na(PA)] <- NA

ndvimax_PA <- ndvimax
ndvimax_PA[is.na(PA)] <- NA

ndvimin_PA <- ndvimin
ndvimin_PA[is.na(PA)] <- NA

pre_PA <- pre
pre_PA[is.na(PA)] <- NA

preD_PA <- preD
preD_PA[is.na(PA)] <- NA

slope_PA <- slope
slope_PA[is.na(PA)] <- NA

soilPH_PA <- soilPH
soilPH_PA[is.na(PA)] <- NA

treeH_PA <- treeH
treeH_PA[is.na(PA)] <- NA

tmax_PA <- tmax
tmax_PA[is.na(PA)] <- NA

treeD_PA <- treeD
treeD_PA[is.na(PA)] <- NA

tseas_PA <- tseas
tseas_PA[is.na(PA)] <- NA



lm_arid <- lm(values(all_PA) ~ values(arid_PA))
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))


summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)


## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.0096396,
               0.0105628,
               0.0248894,
               0.0104246,
               0.0208182,
               0.0197987,
               0.0623993,
               0.0168243,
               0.0140457,
               0.0383225,
               0.0204521)
pente <- c(1,
           1.0198091, 
           0.9962408, 
           0.9611624, 
           1.0133317, 
           0.9854836, 
           0.9609167, 
           0.8840801, 
           0.9764035,
           0.9964867,
           0.9293866,
           0.9638219)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(20:30)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))


##

## Realm 04 - IP

realm <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/realm/Cox_Olson_7realms/Cox_Olson_realm_moll.tif')

arid_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-arid/555_422.tif')
ndvimax_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimax/555_422.tif')
ndvimin_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-ndvimin/555_422.tif')
pre_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-pre/555_422.tif')
preD_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-preD_M/555_422.tif')
slope_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-slope/555_422.tif')
soilPH_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-soilPH/555_422.tif')
treeH_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeH/555_422.tif')
tmax_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tmax_warm_M/555_422.tif')
treeD_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-treeD/555_422.tif')
tseas_PA <- raster('E:/leberro/My Documents/PhD_Paper_1_globalForest/Results/HS_realm/HS_10var/HS_10var-tseas/555_422.tif')

all_PA <- all
all_PA[arid_PA==0] <- NA
plot(all_PA)

lm_arid <- lm(values(all_PA) ~ values(arid_PA))
lm_ndvimax <- lm(values(all_PA) ~ values(ndvimax_PA))
lm_ndvimin <- lm(values(all_PA) ~ values(ndvimin_PA))
lm_pre <- lm(values(all_PA) ~ values(pre_PA))
lm_preD <- lm(values(all_PA) ~ values(preD_PA))
lm_slope <- lm(values(all_PA) ~ values(slope_PA))
lm_soilPH <- lm(values(all_PA) ~ values(soilPH_PA))
lm_treeH <- lm(values(all_PA) ~ values(treeH_PA))
lm_tmax <- lm(values(all_PA) ~ values(tmax_PA))
lm_treeD <- lm(values(all_PA) ~ values(treeD_PA))
lm_tseas<- lm(values(all_PA) ~ values(tseas_PA))


summary(lm_arid) # orange
summary(lm_ndvimax) # darkgreen
summary(lm_ndvimin) # lightgreen
summary(lm_pre)
summary(lm_preD)
summary(lm_slope)
summary(lm_soilPH)
summary(lm_treeH)
summary(lm_tmax)
summary(lm_treeD)
summary(lm_tseas)

## ggplot 2

library(ggplot2)

var <- c('arid','ndvimax','ndvimin','pre','preD_M','slope','soilPH','treeH','tmax_warm','treeD','tseas')

intercept <- c(0,
               0.101061,
               0.112138,
               0.082720,
               0.089272,
               0.090723,
               0.101231,
               0.100255,
               0.103174,
               0.123084,
               0.098577,
               0.090606)
pente <- c(1,
           0.712338, 
           0.691183, 
           0.721455, 
           0.727396, 
           0.716132, 
           0.691556, 
           0.693322, 
           0.684514,
           0.657553,
           0.704424,
           0.715484)

tab1 <- cbind(var,as.numeric(intercept))
tab2 <- cbind(tab1,as.numeric(pente))
tab.df <- as.data.frame(tab2)

head(tab.df)

HS_10var <- c(15:20)
HS_11var <- c(20:30)

xy <- cbind(HS_11var,HS_10var)
xy.df <- as.data.frame(xy)

p <- ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank()
p0 <- p + geom_abline(intercept=intercept[1],slope=pente[1], col='red', show.legend=T, size=1.2, linetype=2) + scale_colour_manual(values="red")
p1 <- p0 + geom_abline(intercept=intercept[2],slope=pente[2], col='orange', show.legend=T, size=0.8, linetype=1) # aridity => 1
p2 <- p1 + geom_abline(intercept=intercept[3],slope=pente[3], col='darkgreen', show.legend=T, size=0.8) # ndvimax => 4
p3 <- p2 + geom_abline(intercept=intercept[4],slope=pente[4], col='lightgreen', show.legend=T, size=0.8) # ndvimin => 9
p4 <- p3 + geom_abline(intercept=intercept[5],slope=pente[5], col='blue', show.legend=T, size=0.8) # pre => 2
p5 <- p4 + geom_abline(intercept=intercept[6],slope=pente[6], col='lightblue', show.legend=T, size=0.8) # preD => 7
p6 <- p5 + geom_abline(intercept=intercept[7],slope=pente[7], col='grey', show.legend=T, size=0.8) # slope => 10
p7 <- p6 + geom_abline(intercept=intercept[8],slope=pente[8], col='brown', show.legend=T, size=0.8) # soilpH => 11
p8 <- p7 + geom_abline(intercept=intercept[9],slope=pente[9], col='turquoise', show.legend=T, size=0.8) # treeH => 6
p9 <- p8 + geom_abline(intercept=intercept[10],slope=pente[10], col='gold', show.legend=T, size=0.8) # tmax => 3
p10 <- p9 + geom_abline(intercept=intercept[11],slope=pente[11], col='purple', show.legend=T, size=0.8) # treeD => 8
p11 <- p10 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8) # tseas => 5

p11

p4

p0 + geom_abline(intercept=intercept[12],slope=pente[12], col='pink', show.legend=T, size=0.8)

# add legend?

p11 + scale_colour_manual(values=c("red",'orange','darkgreen', "blue",'lightblue','grey','brown','turquoise','gold','purple','pink'))

p11 + guides(fill = guide_legend(title = "LEFT", title.position = "left"))

ggplot(data=xy.df,aes(HS_11var,HS_10var)) + geom_blank() + 
  geom_abline(intercept=intercept[1],slope=pente[1], col='red', size=1.2, linetype=2, show.legend=TRUE) + 
  geom_abline(intercept=intercept[2],slope=pente[2], col='orange', size=0.8, linetype=1, show.legend=TRUE) +
  scale_fill_manual(name='My Lines', values=c("black", "blue"))
