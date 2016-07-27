# Method step by step

## I. Raster selection: 
***script folder:*** *1_HS_raster_selection* 

### 1. Variables :

1.	Transformation :

    1. Projection on Mollweide using ArcGis because ESRI extension for Aridity index and Global PET Monthly
      - ‘Project raster’ function from Data Management

    2. Aggregation from 1 to 5 km resolution using GRASS (see section II) 
    
      *inVars_projection.py*
      
    3. Aggregation from 1 to 5 km resolution using R  to have the forest cover for correlation analysis
    
      *variable_res_change.r*

2.	Variables analysis:

    1. Correlation analysis
    
     *variables_cor.R*

      - Variables tested (19):
    
        - BIO1 = Annual Mean Temperature
        - BIO4 = Temperature Seasonality
        - BIO5 = Max Temperature of Warmest Month
        - BIO6 = Min Temperature of Coldest Month
        - BIO7 = Temperature Annual Range (BIO5-BIO6)
        - BIO10 = Mean Temperature of Warmest Quarter
        - BIO11 = Mean Temperature of Coldest Quarter
        - BIO12 = Annual Precipitation
        - BIO14 = Precipitation of Driest Month
        - BIO17 = Precipitation of Driest Quarter
        - Growing degree-days on 0 degree base
        - Aridity index
        - Potential evapotranspiration seasonality
        - NDVI (min) (annual mean) 
        - NDVI (max) (annual mean)
        - Slope
        - Tree density (%) 
        - Canopy height (m)
        - Soil pH for the topsoil - 0 to 5 cm

      -	Variable selection (11):
    
        - tseason: Temperature seasonality
        - tmax_warm_mt: Temperature max of the warmest month
        - pre: Annual precipitation
        - preDry: Precipitation of the driest month
        -	arid: Aridity index
        - ndvimin: NDVI minimum 
        -	ndvimax: NDVI maximum
        -	treeH: Tree height
        -	treeD: Tree cover density
        -	soil: soil pH
        -	slope

2. Realm selection:

  Extract realm ArcGis by using dissolve (merge different polygons in one, e.g.ecoregions).
  
  1.	Mix between Cox and Olson map, using Olson boundaries as Cox digital version doesn’t exist:
    -	Add a column in the attribute of the vector files (edit)
    -	Dissolve the realms

  2.	Creation of 6 realms:
      1.	Neartic (Olson)
      2.	Paleartic (Olson)
      3.	South American (Cox)
      4.	African (Cox)
      5.	Indo-Pacific (Cox) => West delimitation around Himalaya according to forest class
      6.	Australian (Cox)

    Projection in Mollweide using ArcGis (Project raster function from Data Management)

3. Forest class selection: 

*From http://www.esa-landcover-cci.org/*

*See E:/leberro/My Documents/PhD_Paper_1_globalForest/Database/ESA_GLC/ESACCI-LC-Legend.xls*

1-	Tree broadleaved evergreen

2-	Tree broadleaved deciduous

3-	Tree needleleaved evergreen

4-	Tree needleleaved deciduous

5-	Tree mixed leaf type

6-	Tree flooded, fresh water

7-	Tree flooded, saline water

    -	Steps:
        -	Projection in Mollweide using ArcGis (Project raster function from Data Management)

        -	Change resolution at 5km using ArcGis (‘resample’ option using ‘majority ‘option) : 
        *Raster: E: /leberro/My Documents/PhD_Paper_1_globalForest/Database/ESA_GLC/ esa_glc_moll_5km_resample.tif*

        -	Project raster to match the raster with input variables and select forested areas using GRASS rules (r.reclass) : 
        *esa_glc_transformation_reclass.py*
        *Raster: E: /leberro/My Documents/PhD_Paper_1_globalForest/Database/ESA_GLC/ esa_forest_moll5km.tif*

        -	To select 5 km cells which have more that 50 % of forest while aggregation:
            -	Change the extent of esa_glc_moll300m to have a grid which match with our forest raster
                -	By modifying the extent, the resolution is changed as well… So we adapt the extent to the 1km input variable (pet_season)
                -	GGIS: Raster -> Raster calculator, load eas_glc_moll300m and eas_forest_moll5km. Enter esa_glc_moll300m@1, then click on esa_forest_moll5km@1 raster in Raster bands to highlight it and click on Current layer extent button. 
            -	GRASS: Reclassify all 7 forests categories as forest
                Script:
            -	Create a grid based on 5 km resolution raster esa_forest_moll5km.tif (or input variable)
            - 	Overlap the grid of 5km on the raster esa_forest_moll5km.tif 
            -	Select the 5km cells which have more than 50 % tree cover (from class 50 to 90) : aggregate the 0-1 map created further by 5 using sum
            -	Create a raster map of 5km from this selection
            -	Clip it to the current esa_glc_moll5km.tif

II. Preprocessing – HS files preparation

-	inVars: set projection, resolution at 5km, NA values (see gdalinfo rast.tif ) with the script:
script/HS_realm/ 1_HS_raster_selection/inVars_projection.py

-	Pas = ecoregions: each forested realms: Cox_Olson_realm_forest.shp
-	Run pas: 
> python 
>>> from subpas_loop import *

-	Ecoregion = pas: world forests raster: Forest 1 / no Forest 0
o	Copy/paste of the pas output from subpas_loop.py => pa_1 => eco_555
o	To reproject the pas files: 
Script: script/HS_realm/realm_forested_raster.r
-	Run ehab:
> python
>>> from ehab_realm import *
>>> ehabitat(‘555’,’’,’’)

III. Postprocessing – eHab output

postprocessing1_raster_region.py
postprocessing2_merge_rasters.r

To assess the importance of each variable:
HS computed 12 times :
1.	1 for HS_all = computed with all variables (11)
2.	11 for HS_all-1var = computed with 10 variables => 11 variables – 1 variable, the variable we want to assess the importance in the HS computation
3.	 Performing linear regressions with lm( HS_all ~ HS_all-1var)
Script: regression_importance_variables.r
Ex:  
HS_all-aridity => HS was computed with 10 variables, aridity index was excluded
lm (HS_all ~ HS_all-aridity) => by comparing the slope of HS computed with all variables – HS computed with one, assess the importance of the variable missing

IV. Species data processing

Data received from Graeme:
-	12791 files: 
o	esh_spnum_1 : resident birds 
o	esh_spnum_2 : breeding area for migrant
o	esh_spnum_3 : non breeding area for migrant
nb: some sp are also _12, _13, _23, _123

esh_spnum_1 : length(list_1) # 9679 sp => 8366 sp (without _12.tif, _13.tif, _123.tif) => 9679-293-352-668 = 8366
esh_spnum_2  : length(list_2) # 1530 sp => 10 sp (without _23.tif, _123.tif) => 1530-352-668 = 10
esh_spnum_3 : length(list_3) # 1582 sp => 3 sp (without _13.tif, _23.tif, _123.tif) => 1582-352-668 = 3 => error for these species without breeding range
esh_spnum_12  : length(list_12) # 293 sp => species with populations which are partially migrants, sp with both resident range and breeding range for migrants pop.
esh_spnum_13  : length(list_13) # 352 sp 
esh_spnum_23  : length(list_23) # 559 sp 
esh_spnum_123  : length(list_123) # 668 sp 

-	resolution: 5km

-	values from 0 to 1 reflecting probability of occurrence

Data treatment:
Method:
Number of small range species map:
1.	Project _1.tif and _2.tif species in moll, resolution 5 km: ArcGIS using Iterator, see birds_breeding_range.r
 
2.	Change resolution in 100 km : bird_rast_01.py
3.	Select forest species (>= 1/3rd of their breeding range in forest): birds_breeding_range.r
4.	Merge the breeding areas with non-breeding areas of the selected forest species at 5km resolution to have the whole species _1.tif, _2.tif, _3.tif => _1.tif, _2.tif, _12.tif, _13.tif, _123.tif : sp_assemblage.r
5.	Project and change the resolution at 100km for these raster files: ArcGIS using Iterator (see birds_breeding_range.r) and bird_rast_01.py
6.	Select small range size species => first quantile, compare with nb species <= 50 000km2 : birds_breeding_range.r
7.	Change values in 0/1 (binary) of the selected forest species (all ranges, small and big): bird_rast_01.py
8.	Sum the rasters of restricted range species (first quantile and/or 50 000km2)
9.	Sum the rasters of all forest species
10.	Produce the map: [nb restricted range sp]/[forest sp nb]

