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
