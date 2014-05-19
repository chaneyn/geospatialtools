import gdal
import osgeo
import os
from osgeo import osr
import numpy as np

def extract_point_data(file,lats,lons):
 
 #Open file and get geotransformation
 ds = gdal.Open(file)
 gt = ds.GetGeoTransform()
 rb = ds.GetRasterBand(1)

 #Compute ilats and ilons
 ilons = np.round((np.array(lons) - gt[0])/gt[1]).astype(np.int)
 ilats = np.round((np.array(lats) - gt[3])/gt[5]).astype(np.int)

 #Extract data
 values = []
 for i in xrange(ilons.size):
  values.append(rb.ReadAsArray(ilons[i],ilats[i],1,1)[0])

 return np.array(values).astype(np.float)

def read_raster(file):

 #Read in the raster
 dataset = gdal.Open(file)

 #Get dimensons
 nx = dataset.RasterXSize
 ny = dataset.RasterYSize

 #Retrieve band
 band = dataset.GetRasterBand(1)

 #Convert to numpy array
 data = band.ReadAsArray(0,0,nx,ny).astype(np.float32)

 return data

def write_raster(metadata,data,file):

 cols = metadata['nlon']
 rows = metadata['nlat']
 minlon = metadata['minlon']
 if minlon > 180: minlon = minlon - 360
 bands = 1
 driver = gdal.GetDriverByName('GTiff')
 #Create file
 ds = driver.Create(file,cols,rows,1,gdal.GDT_Float32)
 #Set geo information
 ds.SetGeoTransform([minlon,metadata['res'],0,metadata['maxlat'],0,-metadata['res']])
 proj = osr.SpatialReference()
 proj.SetWellKnownGeogCS("EPSG:4326")
 ds.SetProjection(proj.ExportToWkt())
 outband = ds.GetRasterBand(1)
 outband.WriteArray(np.flipud(data),0,0)
 ds = None

 return

def shapefile2raster(raster_in,shp_in,raster_out,workspace,field):

 #Extract coordinates and projection info from the target file
 ds = gdal.Open(raster_in)
 gt = ds.GetGeoTransform()
 cols = ds.RasterXSize
 rows = ds.RasterYSize
 srs = osgeo.osr.SpatialReference()
 srs.ImportFromWkt(ds.GetProjection())
 proj4 = srs.ExportToProj4()

 #Rasterize the shapefile
 shp_out = '%s/%d' % (workspace,np.random.randint(10**5))
 os.system("ogr2ogr -t_srs '%s' %s %s" % (proj4,shp_out,shp_in))
 minx = gt[0]
 miny = gt[3]+rows*gt[5]
 maxx = gt[0]+cols*gt[1]
 maxy = gt[3]
 os.system('gdal_rasterize -init -9999 -a %s -te %.16f %.16f %.16f %.16f -tr %.16f %.16f %s %s' % (field,minx,miny,maxx,maxy,gt[1],gt[5],shp_out,raster_out))
 os.system('rm -r %s' % shp_out)

 return

def raster2raster(raster_template,raster_in,raster_out):

 #Extract coordinates and projection info from the target file
 ds = gdal.Open(raster_template)
 gt = ds.GetGeoTransform()
 cols = ds.RasterXSize
 rows = ds.RasterYSize
 srs = osgeo.osr.SpatialReference()
 srs.ImportFromWkt(ds.GetProjection())
 proj4 = srs.ExportToProj4()

 #Regrid the raster
 minx = gt[0]
 miny = gt[3]+rows*gt[5]
 maxx = gt[0]+cols*gt[1]
 maxy = gt[3]
 os.system("gdalwarp -overwrite -t_srs '%s' --config GDAL_CACHEMAX 5000 -te %.16f %.16f %.16f %.16f -tr %.16f %.16f %s %s" % (proj4,minx,miny,maxx,maxy,gt[1],gt[5],raster_in,raster_out))

 return

