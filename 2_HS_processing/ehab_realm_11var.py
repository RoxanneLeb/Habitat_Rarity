#### Script by Javier Martinez-Lopez (UTF-8)

# import library
from __future__ import division # from a library import function
from datetime import datetime
import numpy as np # as "nickname" we recall the library, np = numpy
import scipy.ndimage as nd
import os.path
import scipy
from scipy.linalg import cholesky, solve_triangular
from scipy.spatial import distance
from scipy.stats import chisqprob
from sklearn.externals.joblib import Parallel, delayed
from multiprocessing import cpu_count
import csv
import os
import sys
from osgeo import ogr,gdal

## fonction which will be used in the script: 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# fonctions to compute Mahalanobis distances in parallel
	# 1. safelyMakeDir : to create new directory
	# 2. _schedule : to control the size of ??
	# 3. _mahalanobis_distances : ?
	# 4. _mahalanobis_distances_scipy: ?
	# 5. initglobalmaps: to import all the input data (variables)
	# 6. ehabitat: ?

# 1. safelyMakeDir	
def safelyMakeDir(d): # Create  directory safely - don't disturb it if it exists already - LB
	try:
		os.makedirs(d) 
		return True
	except OSError:
		if os.path.isdir(d):
			print("Can't create directory: %s - a directory with this name already exists. It will be used for the results of the analysis." % d)	
			return True
		else: 
			print("Can't create directory: %s." % d)	
			return False

# 2. _schedule	# Mahalanobis functions by Sturla Molden

def _schedule(n, nproc):
	"""guided scheduler""" 
	start = 0
	size = (n - start) // nproc
	while size > 100:
		yield slice(start, start + size)
		start += size
		size = (n - start) // nproc
	yield slice(start, n + 1)
	return

# 3. _mahalanobis_distances		# Mahalanobis functions by Sturla Molden
def _mahalanobis_distances(m, L, X):
	cX = X - m[np.newaxis, :]
	tmp = solve_triangular(L, cX.T, lower=True).T
	tmp **= 2
	# return np.sqrt(tmp.sum(axis=1))
	return tmp.sum(axis=1)

def mahalanobis_distances(m, S, X, parallel=True):
	L = cholesky(S, lower=True)
	n = X.shape[0]
	if parallel:
		nproc = cpu_count()
		res = (Parallel(n_jobs= -1)
				(delayed(_mahalanobis_distances)
				  (m, L, X[s, :])
					for s in _schedule(n, nproc)))
		return np.hstack(res)
	else:
		return _mahalanobis_distances(m, L, X)


# 4. _mahalanobis_distances_scipy    # scipy.spatial.distance.mahalanobis for comparison
def _mahalanobis_distances_scipy(m, SI, X):
	n = X.shape[0]
	mahal = np.zeros(n)
	for i in xrange(X.shape[0]):
		x = X[i,:]
		mahal[i] = distance.mahalanobis(x,m,SI)
	return mahal

def mahalanobis_distances_scipy(m, S, X, parallel=True):
	SI = np.linalg.inv(S)
	n = X.shape[0]
	if parallel:
		nproc = cpu_count()
		res = (Parallel(n_jobs=-1)
				(delayed(_mahalanobis_distances_scipy)
				 (m, SI, X[s,:])
				   for s in _schedule(n,nproc)))
		return np.hstack(res)
	else:
		return _mahalanobis_distances_scipy(m, SI, X)

#~~~~~~~~~~~~~~~~~~~~~~~~~		

gmaps = 0
nwpath = ''

def initglobalmaps(): # to load variables at global scale
	
	#	SHARED FOLDER PATH OR LOCAL DIRECTORY
	indir = os.path.join(os.path.sep, nwpath, 'inVars')
	print indir
	tseasf = 'tseas_5km_NAval.tif'
	tmax_warm_Mf = 'tmax_warm_M_5km_NAval.tif' 
	pref = 'pre_5km_NAval.tif'
	preD_Mf = 'preD_M_5km_NAval.tif' 
	aridf = 'arid_5km_NAval.tif'
	ndviminf = 'ndvimin_5km_NAval.tif'
	ndvimaxf = 'ndvimax_5km_NAval.tif'
	treeDf = 'treeD_5km_NAval.tif' 
	treehf = 'treeH_5km_NAval.tif' 
	soilf = 'soilPH_5km_NAval.tif'
	slopef = 'slope_5km_NAval.tif'
	
	#1
	tseasf_globalfile = os.path.join(os.path.sep, indir, tseasf)
	print tseasf_globalfile # reference to the file
	global	src_ds_tseas_global
	src_ds_tseas_global = gdal.Open(tseasf_globalfile) # open the file
	global	tseas_global # to make the variables available to the other functions (in the global environment)
	tseas_global = src_ds_tseas_global.GetRasterBand(1) # import data information (raster values)
	global	gt_tseas_global
	gt_tseas_global = src_ds_tseas_global.GetGeoTransform() # import data information => to get the coordinates of the map (georeference meta date)
	print 'tseas'
	
	#2
	tmax_warm_Mf_globalfile = os.path.join(os.path.sep, indir, tmax_warm_Mf)
	print tmax_warm_Mf_globalfile # reference to the file
	global	src_ds_tmax_warm_M_global
	src_ds_tmax_warm_M_global = gdal.Open(tmax_warm_Mf_globalfile) # open the file
	global	tmax_warm_M_global # to make the variables available to the other functions (in the global environment)
	tmax_warm_M_global = src_ds_tmax_warm_M_global.GetRasterBand(1) # import data information (raster values)
	global	gt_tmax_warm_M_global
	gt_tmax_warm_M_global = src_ds_tmax_warm_M_global.GetGeoTransform() # import data information => to get the coordinates of the map (georeference meta date)
	print 'tmax_warm_M'
	
	#3
	pref_globalfile = os.path.join(os.path.sep, indir, pref)
	global	src_ds_pre_global
	src_ds_pre_global = gdal.Open(pref_globalfile)
	global	pre_global
	pre_global = src_ds_pre_global.GetRasterBand(1)
	global	gt_pre_global
	gt_pre_global = src_ds_pre_global.GetGeoTransform()
	print 'pre'
	
	#4
	preD_Mf_globalfile = os.path.join(os.path.sep, indir, preD_Mf)
	global	src_ds_preD_M_global
	src_ds_preD_M_global = gdal.Open(preD_Mf_globalfile)
	global	preD_M_global
	preD_M_global = src_ds_preD_M_global.GetRasterBand(1)
	global	gt_preD_M_global
	gt_preD_M_global = src_ds_preD_M_global.GetGeoTransform()
	print 'preD_M'
	
	#5
	aridf_globalfile = os.path.join(os.path.sep, indir, aridf)
	global	src_ds_arid_global
	src_ds_arid_global = gdal.Open(aridf_globalfile)
	global	arid_global
	arid_global = src_ds_arid_global.GetRasterBand(1)
	global	gt_arid_global
	gt_arid_global = src_ds_arid_global.GetGeoTransform()
	print 'arid'
	
	#6
	ndviminf_globalfile = os.path.join(os.path.sep, indir, ndviminf)
	global	src_ds_ndvimin_global
	src_ds_ndvimin_global = gdal.Open(ndviminf_globalfile)
	global	ndvimin_global
	ndvimin_global = src_ds_ndvimin_global.GetRasterBand(1)
	global	gt_ndvimin_global
	gt_ndvimin_global = src_ds_ndvimin_global.GetGeoTransform()
	print 'ndvimin'
	
	#7
	ndvimaxf_globalfile = os.path.join(os.path.sep, indir, ndvimaxf)
	global	src_ds_ndvimax_global
	src_ds_ndvimax_global = gdal.Open(ndvimaxf_globalfile)
	global	ndvimax_global
	ndvimax_global = src_ds_ndvimax_global.GetRasterBand(1)
	global	gt_ndvimax_global
	gt_ndvimax_global = src_ds_ndvimax_global.GetGeoTransform()
	print 'ndvimax'
		
	#8
	treehf_globalfile = os.path.join(os.path.sep, indir, treehf)	
	global	src_ds_treeh_global
	src_ds_treeh_global = gdal.Open(treehf_globalfile)
	global	treeh_global
	treeh_global = src_ds_treeh_global.GetRasterBand(1)
	global	gt_treeh_global
	gt_treeh_global = src_ds_treeh_global.GetGeoTransform()
	print 'treeh'
	
	#9
	treeDf_globalfile = os.path.join(os.path.sep, indir, treeDf)
	global	src_ds_treeD_global
	src_ds_treeD_global = gdal.Open(treeDf_globalfile)
	global	treeD_global
	treeD_global = src_ds_treeD_global.GetRasterBand(1)
	global	gt_treeD_global
	gt_treeD_global = src_ds_treeD_global.GetGeoTransform()
	print 'treeD'
		
	#10
	soilf_globalfile = os.path.join(os.path.sep, indir, soilf)
	global	src_ds_soil_global
	src_ds_soil_global = gdal.Open(soilf_globalfile)
	global	soil_global
	soil_global = src_ds_soil_global.GetRasterBand(1)
	global	gt_soil_global
	gt_soil_global = src_ds_soil_global.GetGeoTransform()
	print 'soil'
	
	#11
	slopef_globalfile = os.path.join(os.path.sep, indir, slopef)
	global	src_ds_slope_global
	src_ds_slope_global = gdal.Open(slopef_globalfile)
	global	slope_global
	slope_global = src_ds_slope_global.GetRasterBand(1)
	global	gt_slope_global
	gt_slope_global = src_ds_slope_global.GetGeoTransform()
	print 'slope'	
	print "Global variables imported"
	global	gmaps
	gmaps = 1 # rule to say global maps are loaded

def	ehabitat(ecor,nw,nwpathout): # main fonction (names of ecoregion, network directory input, network directory input) rq: let the network directory input, network directory input empty as we have all the data on the same local folder

	global	nwpath
	if nw=='':
		nwpath = os.getcwd()
	else:
		nwpath = nw
		
		# to create directory :
		#~~~~~~~~~~~~~~~~~~~~
		
	if gmaps == 0: # to check if global maps (input variables) are loaded, if not, do it !
		initglobalmaps()
	if nwpathout=='': 
		#outdir = 'results'	#	ToDo: locally	create	folder "results"	if it	does	not	exist!
		outdir = os.path.join(os.path.sep, os.getcwd(), 'results')
		safelyMakeDir(outdir)
	else:
		#outdir = nwpathout+'/results'	#	SHARED	FOLDER	PATH
		outdir = os.path.join(os.path.sep, nwpathout, 'results')
		safelyMakeDir(outdir)
		#~~~~~~~~~~~~~~~~~~~~~~~
		
	# to create variables	
	tseaspamin = tseaspamax = tmax_warm_Mpamin = tmax_warm_Mpamax = prepamin = prepamax = preD_Mpamin = preD_Mpamax = aridpamin = aridpamax = ndviminpamin = ndviminpamax = ndvimaxpamin = ndvimaxpamax = treehpamin = treehpamax = treeDpamin = treeDpamax = soilpamin = soilpamax = slopepamin = slopepamax = None #=  = treepamin = treepamax = = ndwipamin = ndwipamax
	tseaspamean = tmax_warm_Mpamean = prepamean = preD_Mpamean = aridpamean = ndviminpamean = ndvimaxpamean = treehpamean = treeDpamean = soilpamean = slopepamean = None # = treepamean ndwipamean =
	
	s = nd.generate_binary_structure(2,2)	# pattern to know wgich pixels are agggregated, usefull if we want to work on similar areas	most restrictive pattern for the landscape	patches , for landscape patern analyses
											
	#	LOCAL FOLDER
	csvname1 = os.path.join(os.path.sep, outdir, 'ecoregs_done.csv') # create name of the file
	print csvname1
	if os.path.isfile(csvname1) == False: # if the file ecoregs_done doesn't exist :
		wb = open(csvname1,'a') # create the file
		wb.write('None') # first line
		wb.write('\n') # next line
		wb.close() # close the file
	#	LOCAL FOLDER	
	csvname = os.path.join(os.path.sep, outdir, 'hri_results.csv') # same than before
	print csvname
	if os.path.isfile(csvname) == False:
		wb = open(csvname,'a')
		wb.write('ecoregion wdpaid averpasim hr2aver pxpa hr1insumaver hriaver nfeatsaver lpratio lpratio2 numpszok lpmaxsize aggregation tseaspamin tseaspamax tmax_warm_Mpamin tmax_warm_Mpamax prepamin prepamax preD_Mpamin preD_Mpamax aridpamin aridpamax ndviminpamin ndviminpamax ndvimaxpamin ndvimaxpamax treehpamin treehpamax treeDpamin treeDpamax soilpamin soilpamax slopepamin slopepamax tseaspamean tmax_warm_Mpamean prepamean preD_Mpamean aridpamean ndviminpamean ndvimaxpamean treehpamean treeDpamean soilpamean slopepamean') #      treepamin treepamax treepamean ndwipamin ndwipamax ndwipamean
		wb.write('\n')
		wb.close()
	ef = 'eco_'+str(ecor)+'.tif' # crate a name of a tif file based on the extent of the eco_ (input)
	ecofile = os.path.join(os.path.sep, nwpath, 'ecoregs', ef) # file + path
	#ecofile = os.path.join(os.path.sep, nwpath, os.path.sep,'ecoregs', os.path.sep, ef)
	print ecofile
	avail = os.path.isfile(ecofile) # does this ecoregion file exist in the folder? T/F
	if avail == True: # if the file exist in the path mentionned (if it doesn't exit, it's ignore this ecoregion without error !!): 
		eco_csv = str(ecor)+'.csv' # str = to convert ecor to a string, to have a name
		print eco_csv
		ecoparksf = os.path.join(os.path.sep, nwpath, 'pas', eco_csv)
		#ecoparksf = os.path.join(os.path.sep, nwpath, os.path.sep, 'pas', os.path.sep, eco_csv)
		print ecoparksf
		#ecoparksf = nwpath+'/pas/'+str(ecor)+'.csv'
		src_ds_eco = gdal.Open(ecofile) # for each ecoregion, open the tif
		eco = src_ds_eco.GetRasterBand(1) # band 1 called eco
		eco_mask0 = eco.ReadAsArray(0,0,eco.XSize,eco.YSize).astype(np.int32) # read the values of the raster as a matrix as integers (red-yellow 1/0 raster values ex:555.tif), with the number of cells we want to read 'eco.XSize,eco.YSize' (read total size of rows and columns)
		eco_mask = eco_mask0.flatten() # convert columns and rows of a matrix to one vector (eco.XSize*eco.YSize)
		gt_eco = src_ds_eco.GetGeoTransform() # to get the coordinates, properties of the raster map, to begin by a corner
		print 'eco mask'
		xoff = int((gt_eco[0]-gt_tseas_global[0])/5000)
		yoff = int((gt_tseas_global[3]-gt_eco[3])/5000)
		tseas_eco_bb0 = tseas_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		tseas_eco_bb = tseas_eco_bb0.flatten()
		tseas_eco0 = np.where(eco_mask == 1,(tseas_eco_bb),(0))
		tseas_eco = np.where(tseas_eco0 == 65535.0,	(float('NaN')),(tseas_eco0))
		masktseas = np.isnan(tseas_eco)
		tseas_eco[masktseas] = np.interp(np.flatnonzero(masktseas),	np.flatnonzero(~masktseas),	tseas_eco[~masktseas])
		print 'eco tseas'
		xoff = int((gt_eco[0]-gt_tmax_warm_M_global[0])/5000)
		yoff = int((gt_tmax_warm_M_global[3]-gt_eco[3])/5000)
		tmax_warm_M_eco_bb0 = tmax_warm_M_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		tmax_warm_M_eco_bb = tmax_warm_M_eco_bb0.flatten()
		tmax_warm_M_eco0 = np.where(eco_mask == 1,	(tmax_warm_M_eco_bb),(0))
		tmax_warm_M_eco = np.where(tmax_warm_M_eco0 == 65535.0,	(float('NaN')),(tmax_warm_M_eco0))
		masktmax_warm_M = np.isnan(tmax_warm_M_eco)
		tmax_warm_M_eco[masktmax_warm_M] = np.interp(np.flatnonzero(masktmax_warm_M),	np.flatnonzero(~masktmax_warm_M),	tmax_warm_M_eco[~masktmax_warm_M])
		print 'eco tmax_warm_M'
		xoff = int((gt_eco[0]-gt_pre_global[0])/5000)
		yoff = int((gt_pre_global[3]-gt_eco[3])/5000)
		pre_eco_bb0 = pre_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		pre_eco_bb = pre_eco_bb0.flatten()
		pre_eco0 = np.where(eco_mask == 1,	(pre_eco_bb),(0))
		pre_eco = np.where(pre_eco0 == 65535.0,	(float('NaN')),(pre_eco0))
		maskpre = np.isnan(pre_eco)
		pre_eco[maskpre] = np.interp(np.flatnonzero(maskpre),	np.flatnonzero(~maskpre),	pre_eco[~maskpre])
		print 'eco pre'		
		xoff = int((gt_eco[0]-gt_preD_M_global[0])/5000)
		yoff = int((gt_preD_M_global[3]-gt_eco[3])/5000)
		preD_M_eco_bb0 = preD_M_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		preD_M_eco_bb = preD_M_eco_bb0.flatten()
		preD_M_eco0 = np.where(eco_mask == 1,	(preD_M_eco_bb),(0))
		preD_M_eco = np.where(preD_M_eco0 == 65535.0,	(float('NaN')),(preD_M_eco0))
		maskpreD_M = np.isnan(preD_M_eco)
		preD_M_eco[maskpreD_M] = np.interp(np.flatnonzero(maskpreD_M),	np.flatnonzero(~maskpreD_M),	preD_M_eco[~maskpreD_M])
		print 'eco preD_M'
		xoff = int((gt_eco[0]-gt_arid_global[0])/5000)
		yoff = int((gt_arid_global[3]-gt_eco[3])/5000)
		arid_eco_bb0 = arid_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		arid_eco_bb = arid_eco_bb0.flatten()
		arid_eco0 = np.where(eco_mask == 1,	(arid_eco_bb),(0))
		arid_eco = np.where(arid_eco0 == 65535.0,	(float('NaN')),(arid_eco0))
		maskarid = np.isnan(arid_eco)
		arid_eco[maskarid] = np.interp(np.flatnonzero(maskarid),	np.flatnonzero(~maskarid),	arid_eco[~maskarid])
		print 'eco arid'
		xoff = int((gt_eco[0]-gt_ndvimin_global[0])/5000)
		yoff = int((gt_ndvimin_global[3]-gt_eco[3])/5000)
		ndvimin_eco_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		ndvimin_eco_bb = ndvimin_eco_bb0.flatten()
		ndvimin_eco0 = np.where(eco_mask == 1,	(ndvimin_eco_bb),(0))
		ndvimin_eco = np.where(ndvimin_eco0 == 65535.0,	(float('NaN')),(ndvimin_eco0))
		maskndvimin = np.isnan(ndvimin_eco)
		ndvimin_eco[maskndvimin] = np.interp(np.flatnonzero(maskndvimin),	np.flatnonzero(~maskndvimin),	ndvimin_eco[~maskndvimin])
		print 'eco ndvimin'
		xoff = int((gt_eco[0]-gt_ndvimax_global[0])/5000)
		yoff = int((gt_ndvimax_global[3]-gt_eco[3])/5000)
		ndvimax_eco_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		ndvimax_eco_bb = ndvimax_eco_bb0.flatten()
		ndvimax_eco0 = np.where(eco_mask == 1,	(ndvimax_eco_bb),(0))
		ndvimax_eco = np.where(ndvimax_eco0 == 65535.0,	(float('NaN')),(ndvimax_eco0))
		maskndvimax = np.isnan(ndvimax_eco)
		ndvimax_eco[maskndvimax] = np.interp(np.flatnonzero(maskndvimax),	np.flatnonzero(~maskndvimax),	ndvimax_eco[~maskndvimax])
		print 'eco ndvimax'
		xoff = int((gt_eco[0]-gt_treeh_global[0])/5000)
		yoff = int((gt_treeh_global[3]-gt_eco[3])/5000)
		treeh_eco_bb0 = treeh_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		treeh_eco_bb = treeh_eco_bb0.flatten()
		treeh_eco0 = np.where(eco_mask == 1, (treeh_eco_bb),(0))
		treeh_eco = np.where(treeh_eco0 == 65535.0,	(float('NaN')),(treeh_eco0))
		masktreeh = np.isnan(treeh_eco)
		treeh_eco[masktreeh] = np.interp(np.flatnonzero(masktreeh),	np.flatnonzero(~masktreeh),	treeh_eco[~masktreeh])
		print 'eco treeh'
		xoff = int((gt_eco[0]-gt_treeD_global[0])/5000)
		yoff = int((gt_treeD_global[3]-gt_eco[3])/5000)
		treeD_eco_bb0 = treeD_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		treeD_eco_bb = treeD_eco_bb0.flatten()
		treeD_eco0 = np.where(eco_mask == 1,	(treeD_eco_bb),(0))
		treeD_eco = np.where(treeD_eco0 == 1.8e+308,	(float('NaN')),(treeD_eco0))
		masktreeD = np.isnan(treeD_eco)
		treeD_eco[masktreeD] = np.interp(np.flatnonzero(masktreeD),	np.flatnonzero(~masktreeD),	treeD_eco[~masktreeD])
		print 'eco treeD'
		xoff = int((gt_eco[0]-gt_soil_global[0])/5000) # compute the begining of the ecoregion (x and y): gt_eco = ecoregion map , gt_soil_global = global map, we convert it at the resolution we want
		yoff = int((gt_soil_global[3]-gt_eco[3])/5000)
		soil_eco_bb0 = soil_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32) # as before eco_mask0 but for global, to knoe where and how big is the ecoregion at the global (bounding box = square), astype(np.float32) to keep values
		soil_eco_bb = soil_eco_bb0.flatten() # to put vales in one single vector
		soil_eco0 = np.where(eco_mask == 1,	(soil_eco_bb),(0)) # where you have values of 1 in the ecoregion (= ecoregion present) => take the values of the inVars raster> because bounding box is a square wich contain the ecoregion 
		soil_eco = np.where(soil_eco0 == 65535.0,	(float('NaN')),(soil_eco0)) # where you have 65535 (= NULL from treeD for instence), replace by NaN
		masksoil = np.isnan(soil_eco) # to create new mask, when you have nan => TRUE if not => FALSE
		soil_eco[masksoil] = np.interp(np.flatnonzero(masksoil),	np.flatnonzero(~masksoil),	soil_eco[~masksoil]) # where you have TRUE, you do an interpolation with nearbourhood to get values instead of NaN
		print 'eco soil'
		xoff = int((gt_eco[0]-gt_slope_global[0])/5000)
		yoff = int((gt_slope_global[3]-gt_eco[3])/5000)
		slope_eco_bb0 = slope_global.ReadAsArray(xoff,yoff,eco.XSize,eco.YSize).astype(np.float32)
		slope_eco_bb = slope_eco_bb0.flatten()
		slope_eco0 = np.where(eco_mask == 1,	(slope_eco_bb),(0))
		slope_eco = np.where(slope_eco0 == 65535.0,	(float('NaN')),(slope_eco0))
		maskslope = np.isnan(slope_eco)
		slope_eco[maskslope] = np.interp(np.flatnonzero(maskslope),	np.flatnonzero(~maskslope),	slope_eco[~maskslope])
		print 'eco slope'
		ind_eco0 = np.column_stack((tseas_eco,tmax_warm_M_eco,pre_eco,preD_M_eco,arid_eco,ndvimin_eco,ndvimax_eco,treeh_eco,treeD_eco,soil_eco,slope_eco)) # we make an array putting all the columns of each variable together #    tree_eco,,ndwi_eco
		print 'ecovars stacked'

		print ecoparksf
		pa_list0 = np.genfromtxt(ecoparksf,dtype='int')	# crear este archivo en subpas! read the txt file (eco_csv file = 555.csv), reand pas od one ecoreg
		pa_list = np.unique(pa_list0) # take each pas of the list just once (in case the pa name appears several times)
		n = len(pa_list) # to have the number of pas computed
		for	px in range(0,n): #	0,n

			pa = pa_list[px] # for each pa 
			print pa

			outfile = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'.tif') # prepare the .tif files to be filled
			outfile2 = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'_lp.tif')
			outfile3 = os.path.join(os.path.sep, outdir, str(ecor)+'_'+str(pa)+'_mask.tif')
			#outfile = outdir+'/'+str(ecor)+'_'+str(pa)+'.tif'	#	LOCAL FOLDER
			pa_infile = 'pa_'+str(pa)+'.tif'

			pa4 = os.path.join(os.path.sep, nwpath, 'pas', pa_infile) # pas input 
			#pa4 = os.path.join(os.path.sep, nwpath, os.path.sep, 'pas', os.path.sep, pa_infile)
			print pa4
			#pa4 = nwpath+'/pas/pa_'+str(pa)+'.tif'

			dropcols = np.arange(10,dtype=int) # create vector from 0 to 10 for 11 variables
			done = os.path.isfile(outfile) # outfile is already created
			avail2 = os.path.isfile(pa4) # check if input is available
			if done == False and avail2 == True: # if the files pa4 exists but the HRI computing is not done yet (* see ptt)
				pafile=pa4
				src_ds_pa = gdal.Open(pafile) # open the pa with gdal
				par = src_ds_pa.GetRasterBand(1) # same than previously with ecoreg
				pa_mask0 = par.ReadAsArray(0,0,par.XSize,par.YSize).astype(np.int32)
				pa_mask = pa_mask0.flatten()
				ind = pa_mask >	0 # create index T/F, T if >0.
				go = 1
				sum_pa_mask = sum(pa_mask[ind])# create the sum of the 1 values to know the size of the pa
				if sum_pa_mask < 2: go = 0	#	not	processing	areas	smaller	than	2	pixels, if size too small, don't go
				print sum_pa_mask # print the size of the pa
				sum_pa_mask_inv = len(pa_mask[pa_mask == 0]) # what is the size of 0 pixels inside the pa
				print sum_pa_mask_inv # print size of 0
				print len(pa_mask)
				ratiogeom = 5000 # when bounding box is bigger than pa (due to wwpa database error)
				if sum_pa_mask > 0: ratiogeom = sum_pa_mask_inv/sum_pa_mask # total pixel outside much bigger than pixel in pa ?
				#print ratiogeom
				gt_pa = src_ds_pa.GetGeoTransform() # get coordinates of pa
				xoff = int((gt_pa[0]-gt_pre_global[0])/5000) # compute from the global map xoff, yoff of pa
				yoff = int((gt_pre_global[3]-gt_pa[3])/5000)
				if xoff>=0 and yoff>=0 and go == 1: # xoff>0 to be sure we are not in a boarder of ecoreg and our pa is > 2
					num_bands=src_ds_eco.RasterCount # to have the number of bands in one ecoreg, determined as 1 in the script
					driver = gdal.GetDriverByName("GTiff") # create a new tif file where to store the output (hri)
					dst_options = ['COMPRESS=LZW'] # compression method of the files as many will be processed
					dst_ds = driver.Create(	outfile,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Float32,dst_options) # to create the tif file empty
					dst_ds.SetGeoTransform(	src_ds_eco.GetGeoTransform())					
					dst_ds.SetProjection(src_ds_eco.GetProjectionRef()) # determine the file projection. preparation of the tif file over 
					
					# 1. tseas
					xoff = int((gt_pa[0]-gt_tseas_global[0])/5000)
					yoff = int((gt_tseas_global[3]-gt_pa[3])/5000)
					tseas_pa_bb0 = tseas_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					tseas_pa_bb = tseas_pa_bb0.flatten()
					tseas_pa0 = tseas_pa_bb[ind]
					tseas_pa = np.where(tseas_pa0 == 65535.0,	(float('NaN')),(tseas_pa0))
					mask2tseas = np.isnan(tseas_pa)
					if mask2tseas.all() == True: # if all the variable values in a pa are only NA, put -8 in dropcol
						dropcols[0] = -0
					else: # if there are also non NA values
						tseas_pa[mask2tseas] = np.interp(np.flatnonzero(mask2tseas),	np.flatnonzero(~mask2tseas),	tseas_pa[~mask2tseas])
						tseas_pa = np.random.random_sample(len(tseas_pa),)/5000 + tseas_pa
						print 'pa tseas'

						tseaspamin = round(tseas_pa.min(),2)
						tseaspamax = round(tseas_pa.max(),2)
						tseaspamean = round(np.mean(tseas_pa),2)
						print tseaspamin
						print tseaspamax
						tseasdiff = abs(tseas_pa.min()-tseas_pa.max())
						if tseasdiff < 0.001: dropcols[0] = -0
					
					# 2. tmax_warm_M
					xoff = int((gt_pa[0]-gt_tmax_warm_M_global[0])/5000)
					yoff = int((gt_tmax_warm_M_global[3]-gt_pa[3])/5000)
					tmax_warm_M_pa_bb0 = tmax_warm_M_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					tmax_warm_M_pa_bb = tmax_warm_M_pa_bb0.flatten()
					tmax_warm_M_pa0 = tmax_warm_M_pa_bb[ind]
					tmax_warm_M_pa = np.where(tmax_warm_M_pa0 == 65535.0,	(float('NaN')),(tmax_warm_M_pa0))
					mask2tmax_warm_M = np.isnan(tmax_warm_M_pa)
					if mask2tmax_warm_M.all() == True:
						dropcols[1] = -1
					else:
						tmax_warm_M_pa[mask2tmax_warm_M] = np.interp(np.flatnonzero(mask2tmax_warm_M),	np.flatnonzero(~mask2tmax_warm_M),	tmax_warm_M_pa[~mask2tmax_warm_M])
						tmax_warm_M_pa = np.random.random_sample(len(tmax_warm_M_pa),)/5000 + tmax_warm_M_pa
						print 'pa tmax_warm_M'

						tmax_warm_Mpamin = round(tmax_warm_M_pa.min(),2)
						tmax_warm_Mpamax = round(tmax_warm_M_pa.max(),2)
						tmax_warm_Mpamean = round(np.mean(tmax_warm_M_pa),2)
						print tmax_warm_Mpamin
						print tmax_warm_Mpamax
						tmax_warm_Mdiff = abs(tmax_warm_M_pa.min()-tmax_warm_M_pa.max())
						if tmax_warm_Mdiff < 0.001: dropcols[1] = -1			
						
					# 3. pre
					xoff = int((gt_pa[0]-gt_pre_global[0])/5000)
					yoff = int((gt_pre_global[3]-gt_pa[3])/5000)
					pre_pa_bb0 = pre_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					pre_pa_bb = pre_pa_bb0.flatten()
					pre_pa0 = pre_pa_bb[ind]
					pre_pa = np.where(pre_pa0 == 65535.0,	(float('NaN')),(pre_pa0))
					mask2pre = np.isnan(pre_pa)
					if mask2pre.all() == True:
						dropcols[2] = -2
					else:
						pre_pa[mask2pre] = np.interp(np.flatnonzero(mask2pre),	np.flatnonzero(~mask2pre),	pre_pa[~mask2pre])
						pre_pa = np.random.random_sample(len(pre_pa),)/5000 + pre_pa
						print 'pa pre'

						prepamin = round(pre_pa.min(),2)
						prepamax = round(pre_pa.max(),2)
						prepamean = round(np.mean(pre_pa),2)
						print prepamin
						print prepamax
						prediff = abs(pre_pa.min()-pre_pa.max())
						if prediff < 0.001: dropcols[2] = -2
						
					# # 4. preD_M
					xoff = int((gt_pa[0]-gt_preD_M_global[0])/5000)
					yoff = int((gt_preD_M_global[3]-gt_pa[3])/5000)
					preD_M_pa_bb0 = preD_M_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					preD_M_pa_bb = preD_M_pa_bb0.flatten()
					preD_M_pa0 = preD_M_pa_bb[ind]
					preD_M_pa = np.where(preD_M_pa0 == 65535.0,	(float('NaN')),(preD_M_pa0))
					mask2preD_M = np.isnan(preD_M_pa)
					if mask2preD_M.all() == True:
						dropcols[3] = -3
					else:
						preD_M_pa[mask2preD_M] = np.interp(np.flatnonzero(mask2preD_M),	np.flatnonzero(~mask2preD_M),	preD_M_pa[~mask2preD_M])
						preD_M_pa = np.random.random_sample(len(preD_M_pa),)/5000 + preD_M_pa
						print 'pa preD_M'

						preD_Mpamin = round(preD_M_pa.min(),2)
						preD_Mpamax = round(preD_M_pa.max(),2)
						preD_Mpamean = round(np.mean(preD_M_pa),2)
						print preD_Mpamin
						print preD_Mpamax
						preD_Mdiff = abs(preD_M_pa.min()-preD_M_pa.max())
						if preD_Mdiff < 0.001: dropcols[3] = -3			
						
					# 5. arid
					xoff = int((gt_pa[0]-gt_arid_global[0])/5000)
					yoff = int((gt_arid_global[3]-gt_pa[3])/5000)
					arid_pa_bb0 = arid_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					arid_pa_bb = arid_pa_bb0.flatten()
					arid_pa0 = arid_pa_bb[ind]
					arid_pa = np.where(arid_pa0 == 65535.0,	(float('NaN')),(arid_pa0))
					mask2arid = np.isnan(arid_pa)
					if mask2arid.all() == True:
						dropcols[4] = -4
					else:
						arid_pa[mask2arid] = np.interp(np.flatnonzero(mask2arid),	np.flatnonzero(~mask2arid),	arid_pa[~mask2arid])
						arid_pa = np.random.random_sample(len(arid_pa),)/5000 + arid_pa
						print 'pa arid'

						aridpamin = round(arid_pa.min(),2)
						aridpamax = round(arid_pa.max(),2)
						aridpamean = round(np.mean(arid_pa),2)
						print aridpamin
						print aridpamax
						ariddiff = abs(arid_pa.min()-arid_pa.max())
						if ariddiff < 0.001: dropcols[4] = -4					
					
					# 6. ndvimin
					xoff = int((gt_pa[0]-gt_ndvimin_global[0])/5000)
					yoff = int((gt_ndvimin_global[3]-gt_pa[3])/5000)
					ndvimin_pa_bb0 = ndvimin_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					ndvimin_pa_bb = ndvimin_pa_bb0.flatten()
					ndvimin_pa0 = ndvimin_pa_bb[ind]
					ndvimin_pa = np.where(ndvimin_pa0 == 65535.0,	(float('NaN')),(ndvimin_pa0))
					mask2ndvimin = np.isnan(ndvimin_pa)
					if mask2ndvimin.all() == True:
						dropcols[5] = -5
					else:
						ndvimin_pa[mask2ndvimin] = np.interp(np.flatnonzero(mask2ndvimin),	np.flatnonzero(~mask2ndvimin),	ndvimin_pa[~mask2ndvimin])
						ndvimin_pa = np.random.random_sample(len(ndvimin_pa),)/5000 + ndvimin_pa
						print 'pa ndvimin'

						ndviminpamin = round(ndvimin_pa.min(),2)
						ndviminpamax = round(ndvimin_pa.max(),2)
						ndviminpamean = round(np.mean(ndvimin_pa),2)
						print ndviminpamin
						print ndviminpamax
						ndvimindiff = abs(ndvimin_pa.min()-ndvimin_pa.max())
						if ndvimindiff < 0.001: dropcols[5] = -5
					
					# 7. ndvimax
					xoff = int((gt_pa[0]-gt_ndvimax_global[0])/5000)
					yoff = int((gt_ndvimax_global[3]-gt_pa[3])/5000)
					ndvimax_pa_bb0 = ndvimax_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					ndvimax_pa_bb = ndvimax_pa_bb0.flatten()
					ndvimax_pa0 = ndvimax_pa_bb[ind]
					ndvimax_pa = np.where(ndvimax_pa0 == 65535.0,	(float('NaN')),(ndvimax_pa0))
					mask2ndvimax = np.isnan(ndvimax_pa)
					if mask2ndvimax.all() == True:
						dropcols[6] = -6
					else:
						ndvimax_pa[mask2ndvimax] = np.interp(np.flatnonzero(mask2ndvimax),	np.flatnonzero(~mask2ndvimax),	ndvimax_pa[~mask2ndvimax])
						ndvimax_pa = np.random.random_sample(len(ndvimax_pa),)/5000 + ndvimax_pa
						print 'pa ndvimax'

						ndvimaxpamin = round(ndvimax_pa.min(),2)
						ndvimaxpamax = round(ndvimax_pa.max(),2)
						ndvimaxpamean = round(np.mean(ndvimax_pa),2)
						print ndvimaxpamin
						print ndvimaxpamax
						ndvimaxdiff = abs(ndvimax_pa.min()-ndvimax_pa.max())
						if ndvimaxdiff < 0.001: dropcols[6] = -6
						
					# 8. treeh
					xoff = int((gt_pa[0]-gt_treeh_global[0])/5000) # start to read the trrh cover in the pa
					yoff = int((gt_treeh_global[3]-gt_pa[3])/5000)					
					treeh_pa_bb0 = treeh_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32) 
					treeh_pa_bb = treeh_pa_bb0.flatten()
					treeh_pa0 = treeh_pa_bb[ind] # to read the values inside the pa when we have 1 (higher than 0)
					treeh_pa = np.where(treeh_pa0 == 255.0, (float('NaN')),(treeh_pa0))
					mask2treeh = np.isnan(treeh_pa) # isnan makes an index of T/F, to know which is nan
					if mask2treeh.all() == True: # if all the variable values in a pa are only NA, put -8 in dropcol
						dropcols[7] = -7
					else: # if there are also non NA values
						treeh_pa[mask2treeh] = np.interp(np.flatnonzero(mask2treeh),	np.flatnonzero(~mask2treeh),	treeh_pa[~mask2treeh]) # we do the neirbourood interpolation
						treeh_pa = np.random.random_sample(len(treeh_pa),)/5000 + treeh_pa # add a random noise with an insignificant value to avoid to have all the pixels (from input variables) with the same value and to avoid perfect inverse correlation (between 2 variables ex:treeD/treehs), otherwise Mahalanobis distance can't work
						print 'pa treeh'

						treehpamin = round(treeh_pa.min(),2) #to compute th min of the variable
						treehpamax = round(treeh_pa.max(),2) #to compute the max of the variable
						treehpamean = round(np.mean(treeh_pa),2) #to compute the mean of the variable
						print treehpamin
						print treehpamax
						print treehpamean
						treehdiff = abs(treeh_pa.min()-treeh_pa.max())
						if treehdiff < 0.001: dropcols[7] = -7 # if the difference is too tiny = if the value of the variable don't change => don't use it
					
					# 9. treeD
					xoff = int((gt_pa[0]-gt_treeD_global[0])/5000)
					yoff = int((gt_treeD_global[3]-gt_pa[3])/5000)
					treeD_pa_bb0 = treeD_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					treeD_pa_bb = treeD_pa_bb0.flatten()
					treeD_pa0 = treeD_pa_bb[ind]
					treeD_pa = np.where(treeD_pa0 == 1.8e+308,	(float('NaN')),(treeD_pa0))
					mask2treeD = np.isnan(treeD_pa)
					if mask2treeD.all() == True:
						dropcols[8] = -8
					else:
						treeD_pa[mask2treeD] = np.interp(np.flatnonzero(mask2treeD),	np.flatnonzero(~mask2treeD),	treeD_pa[~mask2treeD])
						treeD_pa = np.random.random_sample(len(treeD_pa),)/5000 + treeD_pa
						print 'pa treeD'

						hpamin = round(treeD_pa.min(),2)
						hpamax = round(treeD_pa.max(),2)
						hpamean = round(np.mean(treeD_pa),2)
						print hpamin
						print hpamax
						hdiff = abs(treeD_pa.min()-treeD_pa.max())
						if hdiff < 0.001: dropcols[8] = -8	
						
					# 10. soil
					xoff = int((gt_pa[0]-gt_soil_global[0])/5000)
					yoff = int((gt_soil_global[3]-gt_pa[3])/5000)
					soil_pa_bb0 = soil_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					soil_pa_bb = soil_pa_bb0.flatten()
					soil_pa0 = soil_pa_bb[ind]
					soil_pa = np.where(soil_pa0 == 65535.0,	(float('NaN')),(soil_pa0))
					mask2soil = np.isnan(soil_pa)
					if mask2soil.all() == True:
						dropcols[9] = -9
					else:
						soil_pa[mask2soil] = np.interp(np.flatnonzero(mask2soil),	np.flatnonzero(~mask2soil),	soil_pa[~mask2soil])
						soil_pa = np.random.random_sample(len(soil_pa),)/5000 + soil_pa
						print 'pa soil'

						soilpamin = round(soil_pa.min(),2)
						soilpamax = round(soil_pa.max(),2)
						soilpamean = round(np.mean(soil_pa),2)
						print soilpamin
						print soilpamax
						soildiff = abs(soil_pa.min()-soil_pa.max())
						if soildiff < 0.001: dropcols[9] = -9
						
					# 11. slope
					xoff = int((gt_pa[0]-gt_slope_global[0])/5000)
					yoff = int((gt_slope_global[3]-gt_pa[3])/5000)
					slope_pa_bb0 = slope_global.ReadAsArray(xoff,yoff,par.XSize,par.YSize).astype(np.float32)
					slope_pa_bb = slope_pa_bb0.flatten()
					slope_pa0 = slope_pa_bb[ind]
					slope_pa = np.where(slope_pa0 == 65535.0,	(float('NaN')),(slope_pa0))
					mask2slope = np.isnan(slope_pa)
					if mask2slope.all() == True:
						dropcols[10] = -10
					else:
						slope_pa[mask2slope] = np.interp(np.flatnonzero(mask2slope),	np.flatnonzero(~mask2slope),	slope_pa[~mask2slope])
						slope_pa = np.random.random_sample(len(slope_pa),)/5000 + slope_pa
						print 'pa slope'

						slopepamin = round(slope_pa.min(),2)
						slopepamax = round(slope_pa.max(),2)
						slopepamean = round(np.mean(slope_pa),2)
						print slopepamin
						print slopepamax
						slopediff = abs(slope_pa.min()-slope_pa.max())
						if slopediff < 0.001: dropcols[10] = -10
	
					cols = dropcols[dropcols>=0] # select the "columns" which are positive ex: if hdiff < 0.001: dropcols[3] = -3 => we not use it
					ind_pa0 = np.column_stack((tseas_pa,tmax_warm_M_pa,pre_pa,preD_M_pa,arid_pa,ndvimin_pa,ndvimax_pa,treeh_pa,treeD_pa,soil_pa,slope_pa)) # stack of all the columns (even the negative)  tree_pa, ndwi_pa,
					ind_pa = ind_pa0[:,cols] # we select from the previoous stack only positive columns for pas
					ind_eco = ind_eco0[:,cols] # we select from the previoous stack only positive columns for ecoreg
					print ind_pa.shape
					hr1sum = hr1insum = indokpsz = pszok = sumpszok = lpratio2 = numpszok = hr1averpa = hr3aver = hr2aver = pszmax = num_featuresaver = lpratio = hr1medianpa = hr1insumaver = pxpa = aggregation = None
					print "PA masked"
					#print ind_pa
					if ind_pa.shape[0]>2 and ind_pa.shape[1]>1: #if we have at least 3 pixels per pa and 2 variables, then we start mahalanobis computation
						Ymean = np.mean(ind_pa,axis=0) #mean of each of the positive variables for one pa
						print 'Max. mean value is '+ str(Ymean.max()) # print the maximum of the mean of the 6 variables
						print "Ymean ok"
						Ycov = np.cov(ind_pa,rowvar=False) # do the covariance among the positive variables
						print 'Max. cov value is '+ str(Ycov.max())
						print "Ycov	ok"
						#mh = mahalanobis_distances(Ymean,	Ycov,	ind_eco,	parallel=False)
						#mh = mahalanobis_distances(Ymean,	Ycov,	ind_eco,	parallel=True)
						mh2 = mahalanobis_distances_scipy(Ymean,	Ycov,	ind_eco,	parallel=True) # to compute the mahalanobis distance (in parallel)
						#mh2 = mahalanobis_distances_scipy(Ymean,	Ycov,	ind_eco,	parallel=False)
						maxmh=mh2.max() # max to check there is no NA
						print 'Max. mh value is '+ str(maxmh)
						print 'Max. mh value is nan: '+ str(np.isnan(maxmh))
						mh = mh2*mh2 #multiply the mahalanobis distance by itself, to make chi2 (because sqrt is needed by mahalanobis)
						print "mh ok" # transform with the chi2 the values in 0/1
						pmh = chisqprob(mh,len(cols)).reshape((eco.YSize,eco.XSize)) # chi2 with mh with 6 variables ( or change chisqprob(mh,6) by: chisqprob(mh,len(cols)) ), and transform the vector in 2D matrix
						# pmhh = np.where(pmh	<=	0.001,None,	pmh) # if the value in pmh is really low, put NA (chi2 values goes from 0 to 1)
						# print "pmh ok"	#	quitar	valores	muy	bajos!
						# pmhhmax = pmhh.max() # max should be close to 1
						# print 'Max. similarity value is '+ str(pmhhmax) 
						# dst_ds.GetRasterBand(1).WriteArray(pmhh) #put values of pmhh in ecoregion map (dst_ds, l.399)
						# dst_ds = None # close the file and save it (with hri inside)
						# hr11 = np.where(pmhh>0,1,0) # 0.5						
						print "pmh ok"	#	quitar	valores	muy	bajos!
						pmhhmax = pmh.max() # max should be close to 1
						print 'Max. similarity value is '+ str(pmhhmax) 
						dst_ds.GetRasterBand(1).WriteArray(pmh) #put values of pmh in ecoregion map (dst_ds, l.399)
						dst_ds = None # close the file and save it (with hri inside)
						hr11 = np.where(pmh>0,1,0) # 0.5
						hr1 = hr11.flatten()
						hr1sum = sum(hr1)
						print 'Number of pixels with similarity higher than 0 is '+str(hr1sum)
						hr1insumaver = hr1insum = 0
						hr1sumaver = hr1sum
						src_ds_sim = gdal.Open(outfile)
						sim = src_ds_sim.GetRasterBand(1)
						gt_sim = src_ds_sim.GetGeoTransform()
						xoff = int((gt_pa[0]-gt_sim[0])/5000)
						yoff = int((gt_sim[3]-gt_pa[3])/5000)
						xextentpa = xoff + par.XSize
						yextentpa = yoff + par.YSize
						xless = sim.XSize - xextentpa
						yless = sim.YSize - yextentpa
						xsize = par.XSize
						ysize = par.YSize
						if xoff>0 and yoff>0 and pmhhmax>0.01 and hr1sum>1 and maxmh!=float('NaN'):#and ratiogeom < 100: #	also	checks	if results	are	not	empty

							# reading the similarity ecoregion without the PA (tmp mask) : for landscape metrics
							os.system('gdal_merge.py '+str(ecofile)+' '+str(pa4)+' -o '+str(outfile3)+' -ot Int32') # gdal tools is needed !!! otherwise => error
							hri_pa_bb03 = sim.ReadAsArray().astype(np.float32)
							hri_pa_bb3 = hri_pa_bb03.flatten()
							
							src_ds_sim2 = gdal.Open(outfile3)
							sim2 = src_ds_sim2.GetRasterBand(1)
							gt_sim2 = src_ds_sim2.GetGeoTransform()
							hri_pa_bb02 = sim2.ReadAsArray().astype(np.int32)
							#hri_pa_bb2 = hri_pa_bb02.flatten()
							hri_pa_bb02_max = hri_pa_bb02.max()
							print 'PA: '+str(pa)
							print 'PA (= max) value from mask = '+str(hri_pa_bb02_max)
							if hri_pa_bb02.shape == hri_pa_bb03.shape:
							 hri_pa02 = np.where(hri_pa_bb02 == pa,0,hri_pa_bb03) # hri_pa_bb02_max


							 if xless < 0: xsize = xsize + xless
							 if yless < 0: ysize = ysize + yless
							 hri_pa_bb0 = sim.ReadAsArray(xoff,yoff,xsize,ysize).astype(np.float32)
							 hri_pa_bb = hri_pa_bb0.flatten()
							 indd = hri_pa_bb > 0
							 hri_pa0 = hri_pa_bb[indd]
							 print 'Total number of pixels with similarity values in PA: '+str(len(hri_pa0))
							 hr1averpa = round(np.mean(hri_pa0[~np.isnan(hri_pa0)]),2) # compute the mean similarity inside a pa
							 #print hr1averpa
							 #hr1medianpa = np.median(hri_pa0[~np.isnan(hri_pa0)])
							 print 'mean similarity in the park is '+str(hr1averpa)
							 #hr1insum = sum(np.where(hri_pa0 >= 0.5,	1,0))	#	use	hr1averpa	as	threshold	instead!						
							 ##hr1inaver = np.where(hri_pa0 >= hr1averpa,	1,0)
							 ##hr1insumaver = sum(hr1inaver)
							 #print hr1insum
							 ##labeled_arrayin, num_featuresin = nd.label(hr1inaver,	structure=s)
							 hr1averr = np.where(hri_pa02 >= hr1averpa,	1,0) # pmhh
							 hr1aver = hr1averr.flatten()
							 print 'Total number of pixels with similarity values in ECO: '+str(sum(hr1aver))
							 labeled_arrayaver, num_featuresaver = nd.label(hr1averr,	structure=s)
							 print 'Nr of similar patches found: '+str(num_featuresaver)
							 if num_featuresaver > 0:
							  lbls = np.arange(1, num_featuresaver+1)
							  psizes = nd.labeled_comprehension(labeled_arrayaver, labeled_arrayaver, lbls, np.count_nonzero, float, 0) #-1
							  pszmax = psizes.max()#-hr1insumaver
							  dst_ds2 = driver.Create(outfile2,src_ds_eco.RasterXSize,src_ds_eco.RasterYSize,num_bands,gdal.GDT_Int32,dst_options)
							  dst_ds2.SetGeoTransform(src_ds_eco.GetGeoTransform())
							  dst_ds2.SetProjection(src_ds_eco.GetProjectionRef())
							  dst_ds2.GetRasterBand(1).WriteArray(labeled_arrayaver)
							  dst_ds2 = None
							  #num_feats = num_features - num_featuresaver
							  hr1sumaver = sum(hr1aver)
							  hr2aver = hr1sumaver #- hr1insumaver , number of pixel in the all ecoregion which have similarity values higher or equal than the average inside a pa
							  pxpa = ind_pa.shape[0] # size of the pa used for computation
							  indokpsz = psizes >= pxpa
							  pszsok = psizes[indokpsz] # NEW
							  sumpszok = sum(pszsok)
							  lpratio=round(float(pszmax/pxpa),2)
							  lpratio2=round(float(sumpszok/pxpa),2)
							  numpszok = len(pszsok)
							  hr3aver = round(float(hr2aver/pxpa),2)
							  aggregation = round(float(hr2aver/num_featuresaver),2)
						#hr2 = hr1sumaver - hr1insumaver
						#print hr2
						#hr3 = float(hr2/ind_pa.shape[0])
						#print hr3
					wb = open(csvname,'a')
					var = str(ecor)+' '+str(pa)+' '+str(hr1averpa)+' '+str(hr2aver)+' '+str(pxpa)+' '+str(hr3aver)+' '+str(num_featuresaver)+' '+str(lpratio)+' '+str(lpratio2)+' '+str(numpszok)+' '+str(pszmax)+' '+str(aggregation)+' '+str(tseasmin)+' '+str(tseasmax)+' '+str(tmax_warm_Mmin)+' '+str(tmax_warm_Mmax)+' '+str(premin)+' '+str(premax)+' '+str(preD_Mmin)+' '+str(preD_Mmax)+' '+str(aridmin)+' '+str(aridmax)+' '+str(ndviminmin)+' '+str(ndviminmax)+' '+str(ndvimaxmin)+' '+str(ndvimaxmax)+' '+str(treehmin)+' '+str(treehmax)+' '+str(treeDmin)+' '+str(treeDmax)+' '+str(soilmin)+' '+str(soilmax)+' '+str(slopemin)+' '+str(slopemax)+' '+str(tseasmean)+' '+str(tmax_warm_Mmean)+' '+str(premean)+' '+str(preD_Mmean)+' '+str(aridmean)+' '+str(ndviminmean)+' '+str(ndvimaxmean)+' '+str(treehmean)+' '+str(treeDmean)+' '+str(soilmean)+' '+str(slopemean) #	exclude	PA!	 '+str(treemin)+' '+str(treemax)+' '+str(treemean)+' '+str(ndwimin)+' '+str(ndwimax)+' '+str(ndwimean)+' 
					wb.write(var)
					wb.write('\n')
					wb.close()
					print "results exported"
					os.system('rm '+str(outfile3))
		wb = open(csvname1,'a')	# where we write the final results	(LOCAL	FOLDER)
		var = str(ecor)
		wb.write(var)
		wb.write('\n')
		wb.close()	
	print "END ECOREG: " + str(ecor)

def	run_batch(): #call the previous functions
#	if __name__ == '__main__':
	#from	datetime	import	datetime
	#import	numpy	as	np
	eco_list0 = np.genfromtxt(nwpath + 'pas' + os.path.sep + 'ecoregs.csv',dtype='int')	# read the pas in the ecoregs.csv, crear	este	archivo	en	subpas!
	eco_list = np.unique(eco_list0) # take the unique values
	#print eco_list
	m = len(eco_list) # length
	for	pm	in	range(0,m): # loop	#3,m	#	0	without	the	negative	ecoregs!
		ecor = eco_list[pm] # select the first ID
		print ecor
		 # comput ehab with this ecoregion # eco=ecoregion, '' to put a path if needed. ex: ehabitat(eco_555,'','')
	print str(datetime.now())
	print "BATCH END"