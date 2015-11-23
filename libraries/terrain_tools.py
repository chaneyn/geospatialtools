#Import all the functions from the fortran library
import numpy as np
import terrain_tools_fortran as ttf

def compute_basin_delineation_nbasins(dem,mask,res,nbasins):

 channel_threshold = 10**6
 #Calculate the d8 accumulation area and flow direction
 (area,fdir) = ttf.calculate_d8_acc(dem,res)
 area[mask == 0] = 0.0
 #Iterate until the number of basins match the desired (bisection)
 max_threshold = np.max(area) - res**2
 min_threshold = max_threshold/1000
 #Calculate number of basins for the two boundaries
 channels = ttf.calculate_channels(area,channel_threshold,max_threshold,fdir)
 min_basins = ttf.delineate_basins(channels,mask,fdir)
 min_nbasins = np.unique(min_basins)[1::].size
 print min_nbasins
 #Min iteration
 channels = ttf.calculate_channels(area,channel_threshold,min_threshold,fdir)
 max_basins = ttf.delineate_basins(channels,mask,fdir)
 max_nbasins = np.unique(max_basins)[1::].size
 for i in xrange(10):
  #Calculate midpoint
  c = (np.log(max_threshold) + np.log(min_threshold))/2
  #Calculate the number of basins for the given threshold
  channels = ttf.calculate_channels(area,channel_threshold,np.exp(c),fdir)
  basins = ttf.delineate_basins(channels,mask,fdir)
  c_nbasins = np.unique(basins)[1::].size
  print min_nbasins,c_nbasins,max_nbasins
  #Determine if we have found our solution
  if c_nbasins == nbasins:
   return basins
  #Create the new boundaries
  if nbasins < c_nbasins:
   min_threshold = np.exp(c)
   channels = ttf.calculate_channels(area,channel_threshold,min_threshold,fdir)
   max_basins = ttf.delineate_basins(channels,mask,fdir)
   max_nbasins = np.unique(max_basins)[1::].size
  else:
   max_threshold = np.exp(c)
   channels = ttf.calculate_channels(area,channel_threshold,max_threshold,fdir)
   min_basins = ttf.delineate_basins(channels,mask,fdir)
   min_nbasins = np.unique(min_basins)[1::].size

 print "Did not converge. Returning the best"
 return basins

def define_hrus(basins,dem,channels):

 nbins = 10
 #Define the unique basins
 ubasins = np.unique(basins)[1::]
 tmp = np.copy(basins)
 tmp[:] = 0
 #Create the dem bins
 for basin in ubasins:
  smask = basins == basin
  #Bin the elevation data
  (hist,bins) = np.histogram(dem[smask],bins=nbins)
  #Place the data
  for ibin in xrange(nbins):
   smask = (basins == basin) & (dem >= bins[ibin]) & (dem < bins[ibin+1])
   tmp[smask] = np.mean(dem[smask])
 import matplotlib.pyplot as plt
 tmp = np.ma.masked_array(tmp,tmp==0)
 plt.imshow(tmp)
 plt.show()

 return 

def calculate_hillslope_properties(hillslopes,dem,basins,res,latitude,
    longitude,depth2channel):

 nh = np.max(hillslopes)
 (eh,ah,bh,lath,lonh,erange,hid,d2c) = ttf.calculate_hillslope_properties(hillslopes,
                               dem,basins,res,nh,latitude,longitude,depth2channel)
 properties = {'elevation':eh,
               'area':ah,
               'basin':bh,
               'latitude':lath,
               'longitude':lonh,
               'range':erange,
               'id':hid,
	       'd2c':d2c,
              }
 #Add aspect,slope,covergence,ids

 return properties

def create_nd_histogram(hillslopes,covariates):

 #Define the mask
 mask = hillslopes > 0
 
 #Initialize the cluster number
 icluster = -1

 #Initialize the hru map
 hrus = np.empty(covariates[covariates.keys()[0]]['data'].shape).astype(np.int32)
 hrus[:] = -9999

 #Iterate through each hillslope making the hrus
 uh = np.unique(hillslopes)[1::]
 for ih in uh:
  mask = hillslopes == ih

  #Define the data and the bins
  bins,data = [],[]
  for var in covariates:
   bins.append(covariates[var]['nbins'])
   data.append(covariates[var]['data'][mask])
  bins = np.array(bins)
  data = np.array(data).T

  #Create the histogram
  H,edges = np.histogramdd(data,bins=bins)
  H = H/np.sum(H)

  #Create a dictionary of hru info
  clusters = {}
  Hfl = H.flat
  for i in xrange(H.size):
   coords = Hfl.coords
   if H[coords] > 0:
    icluster = icluster + 1
    clusters[icluster] = {'pct':H[coords]}
    clusters[icluster]['bounds'] = {}
    for var in covariates:
     key = covariates.keys().index(var)
     clusters[icluster]['bounds'][var] = [edges[key][coords[key]],edges[key][coords[key]+1]]
   Hfl.next()

  #Map the hru id to the grid
  for cid in clusters.keys():
   for id in covariates.keys():
    if covariates.keys().index(id) == 0: string = "(covariates['%s']['data'] >= clusters[%d]['bounds']['%s'][0]) & (covariates['%s']['data'] <= clusters[%d]['bounds']['%s'][1]) & mask" % (id,cid,id,id,cid,id)
    else: string = string +  " & (covariates['%s']['data'] >= clusters[%d]['bounds']['%s'][0]) & (covariates['%s']['data'] <= clusters[%d]['bounds']['%s'][1]) & mask" % (id,cid,id,id,cid,id)
   idx = eval('np.where(%s)' % string)
   hrus[idx] = cid + 1

 return hrus

def create_hillslope_tiles(hillslopes,depth2channel,nbins):

 #Define the clusters for each hillslope
 clusters = np.copy(hillslopes)
 uh = np.unique(hillslopes)[1::]
 for ih in uh:
  mask = hillslopes == ih
  (hist,bins) = np.histogram(depth2channel[mask],bins=nbins)
  for ibin in xrange(nbins):
   smask = mask & (depth2channel >= bins[ibin]) & (depth2channel <= bins[ibin+1])
   clusters[smask] = ibin+1

 #Define the hrus
 hrus = nbins*(hillslopes-1) + clusters
 tmp = np.copy(hrus)
 #Create mapping to clean up hrus
 uhrus = np.unique(hrus)[1::]
 mapping = {}
 for i in xrange(uhrus.size):
  tmp[hrus == uhrus[i]] = i + 1
 hrus = tmp
 hrus[hrus <= 0] = -9999

 return (clusters,hrus)

def calculate_hru_properties(hillslopes,tiles,channels,res,nhillslopes,hrus,depth2channel,slope):

 nhru = np.unique(hrus)[1::].size
 (wb,wt,l,hru_position,hid,tid,hru,hru_area,hru_dem,hru_slope) = ttf.calculate_hru_properties(hillslopes,tiles,channels,nhru,res,nhillslopes,hrus,depth2channel,slope)
 hru_properties = {'width_bottom':wb,
                   'width_top':wt,
                   'hillslope_length':l,
                   'hillslope_position':hru_position,
                   'hillslope_id':hid,
                   'tile_id':tid,
                   'hru':hru,
                   'area':hru_area,
                   'slope':hru_slope,
                   'depth2channel':hru_dem}

 return hru_properties
                       
def cluster_hillslopes(hp,hillslopes,nclusters):

 import sklearn.cluster
 area = hp['area']
 area = (area - np.min(area))/(np.max(area) - np.min(area))
 lats = hp['latitude']
 lats = (lats - np.min(lats))/(np.max(lats) - np.min(lats))
 lons = hp['longitude']
 lons = (lons - np.min(lons))/(np.max(lons) - np.min(lons))
 dem = hp['elevation']
 dem = (dem - np.min(dem))/(np.max(dem) - np.min(dem))
 d2c = hp['d2c']
 d2c = (d2c - np.min(d2c))/(np.max(d2c) - np.min(d2c))
 #X = np.array([lats,lons]).T
 X = np.array([area,]).T
 model = sklearn.cluster.KMeans(n_clusters=nclusters)
 clusters = model.fit_predict(X)+1
 #Assign the new ids to each hillslpe
 hillslopes_clusters = ttf.assign_clusters_to_hillslopes(hillslopes,clusters)
 #Determine the number of hillslopes per cluster
 uclusters = np.unique(clusters)
 nhillslopes = []
 for cluster in uclusters:
  nhillslopes.append(np.sum(clusters == cluster))
 nhillslopes = np.array(nhillslopes)

 return (hillslopes_clusters,nhillslopes)

