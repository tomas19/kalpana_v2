## first step is to add the github repo to the path, so functions can be imported
#import sys
#sys.path.append(r'/rsstu/users/j/jcdietri/DHS-CRCoE-2016-2020/tacuevas/github/Kalpana')
## in this case is not necessary since the module is in the same folder
from downscaling import meshRepLen2raster

'''
This script creates a grass location importing the DEM for downscaling and also creates
a new DEM with same resolution and extend with the size of the mesh triangles. This step
is key for the downscaling and can be run in advance, since only depends on the mesh.
'''

fort14 = r'/home/kalpana/fort.14' ## path of the fort.14 file
epsgIn = 4326 ## CRS for lat/lon
epsgOut = 6543 ## projected CRS for NC
pathOut = r'/home/kalpana/NC9.shp' ## full path of the output shapefile 
grassVer = 8.3
pathRasFiles = 'r/home/kalpana'
rasterFiles = 'ncDEMs_epsg6543'
## in this case we will use the same downscaling raster bounding box as the subdomain

meshRepLen2raster(fort14, epsgIn, epsgOut, pathOut, grassVer, pathRasFiles, rasterFiles, 
                  subDomain=os.path.join(pathRasFiles, rasterFiles), nameGrassLocation=None, createGrassLocation=True, 
                  createLocMethod='from_raster')
