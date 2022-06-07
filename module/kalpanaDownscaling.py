import shutil
import os
import subprocess
import sys
import time
import geopandas as gpd
from kalpanaExport import nc2shp, mesh2gdf

'''
    All functions using GRASS GRIS has an argument 'pkg' which is the grass package.
    The package can't be imported before runing the 'grassEnvVar' function so it's
    not directly imported in this module and will be imported when the functions
    are executed.
'''
def delFiles(listFiles, typeFiles, pkg):
    ''' Delete files form a grass location
        Parameters
            listFiles: list
                list of strings with the name of files to delete
            typeFiles: stromg
                raster, vector (see grass manual for more options)
                https://grass.osgeo.org/grass80/manuals/g.remove.html
        Returns
            None
    '''
    for f in listFiles:
        pkg.run_command('g.remove', flags = 'fb', quiet = True, 
                      type = typeFiles, name = f)

def rastersToList(pathRasters, rasterFiles):
    ''' Create list with full path of rasters. 
        Parameters
            pathRasters: str
                path of the raster files. All must be in the same folder
            rasterFiles: list or str
                name(s) of the raster file(s).
        Return
            rasterFiles: list
                list with the complete path of each raster
    '''
    #### check type of rasterFiles input
    if type(rasterFiles) == str:
        rasterFiles = [rasterFiles]
    # else:
        # print("Warning: if DEM raster resolutions do not match, the aggregate DEM " \
                # "resolution will match the resolution of the first input raster.")
    rasterFiles = [os.path.join(pathRasters, r) for r in rasterFiles]
    return rasterFiles

def grassEnvVar(grassVer):
    ''' Include grass to the environmental variables
        Parameters
            grassVer: float
                Version of the grass software (The code was writen for v8.0).
        Return
            None
    '''
    ########### Launch GRASS
    grassBin = f'C:\\Program Files\GRASS GIS {grassVer}\\grass{10*grassVer:0.0f}.bat'
    startCmd = [grassBin, '--config', 'path']

    p = subprocess.Popen(startCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()
    out = out.strip().decode('utf-8')

    if p.returncode != 0:
        print(f'ERROR: {err}', file=sys.stderr)
        print(f'ERROR: Cannot find GRASS GIS {grassVer} start script: {startCmd}', file=sys.stderr)
        sys.exit(-1)

    ########### Add environment variables
    gisbase = out.strip(os.linesep)
    os.environ['GISBASE'] = gisbase

    ########### Init GRASS environment
    gpydir = os.path.join(gisbase, "etc", "python")
    sys.path.append(gpydir)

def createGrassLoc(grassVer, locPath, createLocMethod, myepsg, rasFile):
    ''' Create a Grass location for the downscaling methods. 
        Parameters
            grassVer: float
                Version of the grass software (The code was writen for v8.0).
            locPath: str
                path where the location will be created
            createLocMethod: str
                Two options "from_epsg" (default) or "from_raster" otherwise an error will be thrown.
            myepsg: int
                epsg code (https://epsg.io/), only used if createLocMehot == "from_epsg".
            rasFile: str
                complete path of the raster which will be used to define the coordinate system
        Return
            None
    '''
    ########### Create location
    grassBin = f'C:\\"Program Files"\\"GRASS GIS {grassVer}"\\grass{10*grassVer:0.0f}.bat'

    #### check if location already exist
    if os.path.isdir(locPath):
        shutil.rmtree(locPath)
    else:
        pass
        # os.mkdir(locPath)
        # print(f'{locPath} created')

    #### epsg can be obtained from a raster or specified by the user
    if createLocMethod == 'from_epsg':
        startCmd = grassBin + ' -c epsg:' + str(myepsg) + ' -e ' + locPath
    elif createLocMethod == 'from_raster':
        #print(f'Projection from {rasFile}')
        startCmd = grassBin + ' -c ' + rasFile + ' -e ' + locPath
    else:
        sys.exit('The create location method specified is incorrect. See the docstring!')

    p = subprocess.Popen(startCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    if p.returncode != 0:
        print(f'ERROR: {err}', file = sys.stderror)
        print(f'"ERROR: Cannot create location {startCmd}', file = sys.stderror)
        sys.exit(-1)
    #else:
        #print(f'Created location {locPath}')

def initGrass(locPath, pkg, mapset='PERMANENT'):
    ''' Set grass working environment
        Parameters
            locPath: str
                path of the grass location
            mapset: str. Default PERMANENT
                mapset where the core data of the project can be stored
    '''
    #### init grass environment
    pkg.init(f'{locPath}\\{mapset}')

def importRasters(rasFiles, createLocMethod, pkg):
    ''' Import rasters using Grass gis.
        Parameters
            rasFiles: list
                list with complete path of each raster file
            createLocMethod: str
                Two options "from_epsg" (default) or "from_raster" otherwise an error will be thrown.
        Returns
            rasterOutList: list
                list with imported raster files
    '''
    #import grass.script as gs
    rasterOutList = []
    #### using first raster resolution
    if createLocMethod == 'from_raster':
        for ras in rasFiles:
            rasObj = os.path.basename(ras).split('.')[0]

            pkg.run_command('r.import', overwrite = True, input = ras, output = rasObj)
            rasterOutList.append(rasObj)

    elif createLocMethod == 'from_epsg':
        for ras in rasFiles:
            rasObj = os.path.basename(ras).split('.')[0]
            try:
                pkg.run_command('r.in.gdal', overwrite = True, input = ras, output=rasObj, 
                                flags = 'o')
            except:
                ## check this line, the code didn't work when a large number of rasters
                ## were imported
                time.sleep(10)
                pkg.run_command('r.in.gdal', overwrite=True, input=ras, output=rasObj, 
                                flags = 'o')
            rasterOutList.append(rasObj)

    else:
        print('Raster resolution options are "align" or an integer value')
        sys.exit(-1)
        
    return rasterOutList

def setDownscalingDEM(rasList, pkg, lim=500):
    ''' Create grass region based in input rasters
        Parameters
            rasList: list
                list of the raster objects imported
            lim: int. Default 500
                maximum number of raster to be imported at once before patching them
        Return
            res: float
                downscaling dem resolution in map units
    '''
    ## new line
    pkg.run_command('g.region', raster = rasList, quiet = True)
    
    if 1 < len(rasList) < lim:
        pkg.run_command('r.patch', input = rasList, output = 'dem',
                        overwrite = True, quiet = True)
        delFiles(rasList, 'raster', pkg)
    
    elif len(rasList) >= lim:
        pDems = []
        for i in range(int(len(grassRasList)/lim) + 1):
            sublist = grassRasList[i*lim:(i+1)*lim]
            pkg.run_command('r.patch', input = sublist, output = f'dem{i:03d}',
                overwrite = True, quiet = True)
            pDems.append(f'dem{i:03d}')
            delFiles(sublist, 'raster', pkg)
        pkg.run_command('r.patch', input = pDems, output = 'dem',
                overwrite = True, quiet = True)
        delFiles(pDems, 'raster', pkg)
    
    else:
        pkg.run_command('g.rename', raster = (rasList[0], 'dem'), 
                        overwrite = True, quiet = True)

def vertUnitConvert(conv, pkg):
    ''' Convert vertical units from meter 2 feet or feet 2 meter.
        Parameters
            conv: str
                "m2ft" or "ft2m"
        Return
            None
    '''
    # #import grass.script as gs
    demOrg = f"dem_{conv.split('2')[0]}"
    pkg.run_command('g.rename',raster = f"dem, {demOrg}")
    
    if conv == 'm2ft':
        pkg.mapcalc("$output=if(!isnull(demOrg),demOrg*3.2808399,null())",
                        output="dem",
                        overwrite=True,
                        quiet=True)
    
    elif conv == 'ft2m':
        pkg.mapcalc("$output=if(!isnull(demOrg),demOrg/3.2808399,null())",
                        output="dem",
                        overwrite=True,
                        quiet=True)
    else:
        sys.exit('Wrong argument. Available options are: m2ft or ft2m')

def setGrassEnv(grassVer, pathGrassLocation, createGrassLocation, pkg0, pkg1,
                pathRasFiles, rasterFiles, createLocMethod, myepsg):
    ''' Initialize grass environment
        Parameters
            grassVer: float
                Version of the grass software (The code was writen for v8.0).
            pathGrassLocation: str
                path and name of the grass location
            createGrassLocation: boolean
                True for creating a new location and loading DEMs, false to use an existing location with DEMs already imported
            pathRasFiles: str
                folder of the raster files
            rasterfiles: str or list
                name or names of the raster files to import
            createLocMethod: str
                Two options "from_epsg" (default) or "from_raster" otherwise an error will be thrown.
            myepsg: int
                psg code (https://epsg.io/), only used if createLocMehot == "from_epsg".
    '''
    if createGrassLocation == True:
        ta = time.time()
        rasFiles = rastersToList(pathRasFiles, rasterFiles) ## list with path of rasters to import
        print(f'        rasters to list: {(time.time() - ta)/60: 0.3f} min')
        
        ta = time.time()
        createGrassLoc(grassVer, pathGrassLocation, createLocMethod, myepsg, rasFiles[0])
        print(f'        create location: {(time.time() - ta)/60: 0.3f} min')
        
        ta = time.time()
        initGrass(pathGrassLocation, pkg1)
        print(f'        init grass: {(time.time() - ta)/60: 0.3f} min')
        
        ta = time.time()
        grassRasList = importRasters(rasFiles, createLocMethod, pkg0) ## try rasters with different resolution? check if rasters are in different crs that location
        print(f'        import raster: {(time.time() - ta)/60: 0.3f} min')
        
        ta = time.time()
        setDownscalingDEM(grassRasList, pkg0)
        print(f'        set downscaling dem: {(time.time() - ta)/60: 0.3f} min')
    
    else:
        ta = time.time()
        initGrass(pathGrassLocation, pkg1)
        print(f'        init grass: {(time.time() - ta)/60: 0.3f} min')
        
def setupGrowing(kalpanaShp, attrCol, mesh2ras, meshFile, minArea, pkg):
    ''' Preprocess kalpana shape file
        Parameters
            kalpanaShp: str
                path of the shape file with AdCirc results
            attrCol: str
                name of the attribute column
            mesh2ras: boolean
                True for converting the mesh vector layer to raster
            meshFile: str
                path of the shapefile with the mes elements or the rasterized version of it.
                if mesh2ras is False, meshFile must be the path of a raster, if True of a shapefile.
            minArea: int
                Minimum size of area to be imported (square meters) 
                
    '''
    ## import shape file with max water level
    # pkg.run_command('v.import', input = kalpanaShp, overwrite = True, 
                   # quiet = True, snap = 0.000001)
    t0 = time.time()
    pkg.run_command('v.in.ogr', input = kalpanaShp, overwrite = True,
                    quiet = True, snap = 0.000001, min_area = 10,
                    flags = 'o')
    print(f'        Import kalpana shapefile: {(time.time() - t0)/60:0.2f} min')
    
    t0 = time.time()
    pkg.run_command('v.to.rast', input = os.path.basename(kalpanaShp[:-4]), 
                    type = 'area', output = 'kalpanaRast', use = 'attr', 
                    quiet = True, attribute_column = attrCol, overwrite = True)
    print(f'        Kalpana shape to raster: {(time.time() - t0)/60:0.2f} min')
    
    if mesh2ras == True: #exportMesh is True so meshFile is a shapefile
    
        t0 = time.time()
        # meshShp = os.path.join(pathaux, os.path.splitext(meshFile)[0]+'.shp')
        pkg.run_command('v.in.ogr', input = meshFile, overwrite = True,
                        quiet = True, min_area = int(minArea/100)*100, flags = 'o', snap = 0.000001)
        print(f'        Import mesh shapefile: {(time.time() - t0)/60:0.2f} min')
        
        t0 = time.time()
        pkg.run_command('v.to.rast', input = os.path.splitext(os.path.basename(meshFile))[0], 
                        type = 'area', use = 'attr', quiet = True, 
                        attribute_column = 'repLen', overwrite = True,
                        output = os.path.splitext(os.path.basename(meshFile))[0])
        print(f'        Mesh shape to raster: {(time.time() - t0)/60:0.2f} min')
        t0 = time.time()
        pkg.run_command('r.out.gdal', input = os.path.splitext(os.path.basename(meshFile))[0], 
                        flags = 'm', format = 'GTiff', nodata = -9999, overwrite = True,
                        output = os.path.splitext(meshFile)[0] + '.tif')
        print(f'        Mesh exported as raster: {(time.time() - t0)/60:0.2f} min')
        
    else: #exportMesh is True so meshFile is a raster
        importRasters([meshFile], 'from_raster', pkg)
        
    # reg = pkg.read_command('g.region', flags = 'g').split()
    
    # ewRes = float(reg[7].split("=")[1])
    # nsRes = float(reg[6].split("=")[1])
    # pkg.run_command('g.region', raster = 'dem', nsres = nsRes, 
                   # ewres = ewRes, overwrite = True, quiet = True)
    # pkg.run_command('r.mask', raster = 'dem', quiet = True, 
                   # overwrite = True)
                   
def staticGrow(repLenFactor, pkg, meshFile):
    ''' Execture r.grow.distance algorithm
        Parameters
            repLenFactor: float. Default 0.5
                factor of the representative length of each triangle used as a maximum
                growing distance. Adcirc results are expanded until the distance to the nearest non-null cel
                is equals to the repLenFactor times the representative length of the triangle of the grown cell is in.
            meshFile: str
                name of the rasterized layer with the mesh elements representative size.
    '''
    # if growRadius == 0: ## grow raster cells without a limit distance
        # pkg.run_command('r.grow.distance', input = 'kalpanaRast', value = 'grownRastVal', overwrite = True,
                        # quiet = True)
        
        # pkg.mapcalc("$output = if(!isnull($input), $input, if($base < $new, $new, null()))",
                    # output = 'grownKalpanaRast0', input = 'kalpanaRast', base = 'dem', 
                    # new = 'grownRastVal', quiet = True, overwrite = True)
     
    # elif growRadius > 0: ## grow raster cells with a limiting distance on growRadius
    # growRadiusSq = growRadius**2 #distance on r.grow is calculated using squared metric so it is the squared of the actual distance
    
    pkg.run_command('r.grow.distance', input = 'kalpanaRast', metric = 'squared', distance = 'grownRastDist',
                    value = 'grownRastVal', overwrite = True, quiet = True) ## flag m: distance in meters
    
    repLenFactorSq = repLenFactor**2
    pkg.mapcalc("$output = if(!isnull($input), $input, if($dist <= $fac * $radius * $radius, $new, null()))",
          output = 'grownKalpanaRast0', input = 'kalpanaRast', radius = meshFile, base = 'dem', fac = repLenFactorSq,
          new = 'grownRastVal', dist = 'grownRastDist', quiet = True, overwrite = True)
    # pkg.mapcalc("$output = if(!isnull($input), $input, if($dist <= $radius && $base < $new, $new, null()))",
              # output = 'grownKalpanaRast0', input = 'kalpanaRast', radius = growRadiusSq, base = 'dem',
              # new = 'grownRastVal', dist = 'grownRastDist', quiet = True, overwrite = True)
    
    # else: #growradius < 0 --> shrink results, not sure why it is usefull but it was in the original version of kalpana
        # kv = pkg.region()
        # scale = np.sqrt(kv['nsres'] * kv['ewres'])
        # growRadiusPixelsSq = (growRadius * scale)**2
        # pkg.run_command('r.grow.distance', input = 'kalpanaRast', metric = 'squared', distance = 'grownRastDist', 
                          # value = 'grownRastVal', flags = 'n', quiet = True,
                          # overwrite = True)  ## flag n: distance to the nearest NO-NULL cell. value in pixels
        
        # pkg.mapcalc("$output = if($dist < $radius, null(), $input)", output = 'grownKalpanaRast0', 
                    # radius = growRadiusPixelsSq, input = 'kalpanaRast', dist = 'grownRastDist', 
                    # quiet = True, overwrite = True)
                   

def clumping(rasterGrown, rasterOrg, rasterNew, clumpSizeThreshold, pkg):
    ''' Function to deal with clump of disconnected cells generated after remove cells
        where the ground level is larger than the grown adcirc value
        Parameters
            rasterGrown: str
                name of the grown raster to analize
            rasterOrg: str
                name of the original raster with the non grown adcirc results
            rasterNew: str
                name of the outpur raster
            clumpSizeThreshold: int or float
                threshold to delete isolated cells. Groups of cells smaller than threshold
                will be removed if are not connected to the main water area.
    '''
    pkg.mapcalc("$output = if(!isnull($input), -1, null())", output = 'temp1', 
                   input = rasterGrown, quiet = True, overwrite = True)

    pkg.run_command('r.clump', input = 'temp1', output = 'temp2', 
                   quiet = True, overwrite = True)

    areas = pkg.read_command('r.stats', input = 'temp2', sort = 'desc', 
                            flags = 'c', quiet = True).split()
    
    rastInfo = pkg.parse_command('r.info', map = 'dem', flags = 'g', delimiter = '=')
    nsres = float(rastInfo['nsres'])
    ewres = float(rastInfo['ewres'])
    res = 0.5*(nsres + ewres)
    
    clumpThres = int(clumpSizeThreshold / res)
    
    #if clumpSizeThreshold != 'max':
    ncells = [x for x in areas[3::2] if int(x) >= clumpThres]
    clumpsID = [x for x, y in zip(areas[2::2], areas[1::2]) if int(y) >=  clumpThres]
    #Interleave the two lists 
    areas = [val for tup in zip(clumpsID, ncells) for val in tup]

    reclassList = ''.join([f"{i} = -1\n" for i in areas[::2]])
        
    #else:
     #   reclassList = f"{areas[2]} = -1"
    
    pkg.write_command('r.reclass', input = 'temp2', output = 'temp3', rules = '-', 
                stdin = reclassList, quiet = True, overwrite = True)
    
    # Passes back grown ADCIRC cells if they coincide with the assigned value (-1).
    pkg.mapcalc("$output = if($A == -1, $B, null())", output = rasterNew, A = 'temp3', 
               B = rasterGrown, quiet = True, overwrite = True)

def postProcessStatic(compAdcirc2dem, floodDepth, kalpanaShp, clumpThreshold, pkg, ras2vec=False):
    ''' Postprocess the grown raster
        Parameters
            compAdcirc2dem: boolean
                True for removing ADCIRC cells with values lower than the dem
            flooDepth: boolean
                True for transform water levels to water depth. False for export
                water levelss
            kalpanaShp: str
                path of the shape file with AdCirc results. It is used to save the
                downscaled shp and raster files
            clumpThreshold: int or float
                threshold to delete isolated cells. Groups of cells smaller than threshold will be removed 
                if are not connected to the main water area.
            ras2vec: boolean. Default False
                If False output raster are not saved as shapefiles
    '''

    # if compArcirc2dem is True means ADCIRC cells with values lower than the dem
    # elevation will be removed
    if compAdcirc2dem == True:
        ta = time.time()
        pkg.mapcalc("$output = if(!isnull($dem), if($dem > $adcirc, null(), $adcirc), $adcirc)",
                    output = 'grownKalpanaRast1', adcirc = 'grownKalpanaRast0', 
                    dem = 'dem', quiet = True, overwrite = True)
        print(f'        Delete ground level: {(time.time() - ta)/60:0.3f}')
    else:
        ta = time.time()
        pkg.run_command('g.rename', raster = ('grownKalpanaRast0', 'grownKalpanaRast1'), 
                      overwrite = True, quiet = True)
        print(f'        Rename: {(time.time() - ta)/60:0.3f}')
                      
    ta = time.time()
    clumping('grownKalpanaRast1', 'kalpanaRast', 'grownKalpanaRastLevel', clumpThreshold, pkg)
    print(f'        Delete unconnected clumps: {(time.time() - ta)/60:0.3f}')
    
    
    pathOut = os.path.dirname(kalpanaShp)
    fileOut = os.path.basename(kalpanaShp)[:-4]# remove extension
    
    ta = time.time()
    # export raster as tif
    pkg.run_command('r.out.gdal', input = 'grownKalpanaRastLevel', flags = 'm', format = 'GTiff', nodata = -9999, 
                   output = os.path.join(pathOut, f'{fileOut}_level_downscaled.tif'), 
                   overwrite = True)
    print(f'        export as tif level: {(time.time() - ta)/60:0.3f}')
    
    if ras2vec == True:
        ta = time.time()
        #Export to ESRI shapefile
        pkg.run_command('r.to.vect', input = 'grownKalpanaRastLevel', output = 'grownKalpanaVectLevel', type = 'area')
        print(f'        ras to vector: {(time.time() - ta)/60:0.3f}')
        
        ta = time.time()
        pkg.run_command('v.out.ogr', input = 'grownKalpanaVectLevel', 
                          output = os.path.join(pathOut, f'{fileOut}_level_downscaled'), 
                          type = 'area', format = 'ESRI_Shapefile', flags = 'se', quiet = True, overwrite = True)
        print(f'        export as shp level: {(time.time() - ta)/60:0.3f}')
    
    if floodDepth == True: ## export water depth in flooded cells
        ta = time.time()
        pkg.mapcalc("$output = if(!isnull($dem) && !isnull($grown), $grown - $dem, null())", grown = 'grownKalpanaRastLevel', 
                    dem = 'dem', output = 'grownKalpanaRastDepth', quiet = True, overwrite = True)
        print(f'        compute water depth: {(time.time() - ta)/60:0.3f}')
        
        ta = time.time()
        # export raster as tif
        pkg.run_command('r.out.gdal', input = 'grownKalpanaRastDepth', flags = 'm', format = 'GTiff', nodata = -9999, 
                       output = os.path.join(pathOut, f'{fileOut}_depth_downscaled.tif'), 
                       overwrite = True)
        print(f'        export as tif depth: {(time.time() - ta)/60:0.3f}')
        
        ta = time.time()
        if ras2vec == True:
            #Export to ESRI shapefile
            pkg.run_command('r.to.vect', input = 'grownKalpanaRastDepth', output = 'grownKalpanaVectDepth', type = 'area')
            print(f'        ras to shp depth: {(time.time() - ta)/60:0.3f}')
            
            ta = time.time()
            pkg.run_command('v.out.ogr', input = 'grownKalpanaVectDepth', 
                            output = os.path.join(pathOut, f'{fileOut}_depth_downscaled'), 
                            type = 'area', format = 'ESRI_Shapefile', flags = 'se', quiet = True, overwrite = True)
            print(f'        export as shp depth: {(time.time() - ta)/60:0.3f}')
        
def runStatic(ncFile, levels, epsgIn, epsgOut, vUnitIn, vUnitOut, vDatumIn, vDatumOut, pathOut,  grassVer, pathRasFiles, rasterFiles,
              var='zeta_max', conType ='polygon', subDomain=None, vDatumPath=None, exportMesh=False, n=-1, rs=42, aggfunc='mean',
              pathGrassLocation=None, createGrassLocation=True, createLocMethod='from_raster', attrCol='zMean', repLenGrowing=0.5, 
              meshFile=None, compAdcirc2dem=True, floodDepth=True, clumpThreshold='from_mesh', perMinElemArea=1, ras2vec=False):
    ''' Run static downscaling method and the nc2shp function of the kalpanaExport module.
        Parameters
        ********************************************************************************************************************
        ****************************** REQUIRED inputs of nc2shp function **************************************************
        ********************************************************************************************************************
            ncFile: string
                path of the adcirc output, must be a netcdf file
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            epsgIn: int
                coordinate system of the adcirc input
            epsgOut: int
                coordinate system of the output shapefile
            pathout: string
                complete path of the output file (*.shp)
        ********************************************************************************************************************
        ***************************************** REQUIRED inputs of static method *****************************************
        ********************************************************************************************************************
            grassVer: float
                Version of the grass software (The code was writen for v8.0).
            pathRasters: str
                path of the raster files
            rasterFiles: list or str
                name(s) of the raster file(s).
        ********************************************************************************************************************
        ***************************************** OPTIONAL inputs of nc2shp function ***************************************
        ********************************************************************************************************************
            var: string. DEFAULT zeta_max
                Name of the variable to export
            conType: string. DEFAULT polygon
                'polyline' or 'polygon'
            vDatumPath: string. Default None
                full path of the instalation folder of vdatum (https://vdatum.noaa.gov/). Required only if vertical datums
                vDatumIn and vDatumOut are different
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates. The crs must be the same of the
                adcirc input file.
            exportMesh: boolean. Default False
                True for export the mesh geodataframe and also save it as a shapefile
            n: int
                factor used to compute the number of clusters. n_clusters = n_elements / n
            rs: int. Default 42
                Random state
            aggfunc: str. Default mean
                how to aggregate quantitative values
        ********************************************************************************************************************
        ***************************************** OPTIONAL inputs of static method *****************************************
        ********************************************************************************************************************
            pathGrassLocation: str. DEFAULT None
                path and name of the grass location. If None the grass location will be called 'grassLoc' and save in the 
                same path of the extracted shape file (pathout).
            createGrassLocation: boolean. DEFAULT True
                True for creating a new location and loading DEMs, false to use an existing location with DEMs already imported
             createLocMethod: str. DEFAULT 'from_raster'
                Two options "from_epsg" (default) or "from_raster" otherwise an error will be thrown.
            attrCol: str. DEFAULT 'avgVal'
                name of the attribute column
            repLenFactor: float. Default 0.5
                factor of the representative length of each triangle used as a maximum
                growing distance. Adcirc results are expanded until the distance to the nearest non-null cel
                is equals to the repLenFactor times the representative length of the triangle of the grown cell is in.
            meshFile: string
                path of the raster file with the mesh elements. This file does not depend of the simulation it self, only on the mesh.
                So it is not necessary to export the adcirc mesh as a shapefile each time a simulation is downscaled. For speed up the
                downscaling process this input should be not None only the first time a mesh is used.
            compAdcirc2dem: boolean. DEFAULT True
                True for removing ADCIRC cells with values lower than the dem.
            floodDepth: boolean . DEFAULT True
                True for transform water levels to water depth. False for export water levelss.
            clumpThreshold: str or int or float. Default 'from_mesh'
                threshold to delete isolated cells. Groups of cells smaller than threshold will be removed if are not connected to the main water area. 
                If equals 'from_mesh' the minimum area of the triangles elements is used to define the threshold. If int or float, this value will be used as
                threshold, the unit of this number must be in the map units
            perMinElemArea: float. Default 1
                percentage of the smaller element's area used as a threshold to delete unconnected groups of raster cells. Only used if
                clumpThreshold is "from_mesh".
            ras2vec: boolean. Default False
                For speed up the process is recommended that raster files should not be converted to shapefiles (False).
    '''
    t0 = time.time()
    pathaux = os.path.dirname(pathOut)
    if not os.path.exists(pathaux):
        os.mkdir(pathaux)
    
    if exportMesh == True:
        gdf, mesh = nc2shp(ncFile, var, levels, conType, epsgIn, epsgOut, vUnitIn, vUnitOut, vDatumIn, vDatumOut, pathOut, 
                           vDatumPath, subDomain, exportMesh, n, rs, aggfunc, os.path.splitext(os.path.basename(meshFile))[0])
        meshFile = os.path.join(pathaux, os.path.splitext(os.path.basename(meshFile))[0] + '.shp')
    else:
        gdf = nc2shp(ncFile, var, levels, conType, epsgIn, epsgOut, vUnitIn, vUnitOut, vDatumIn, vDatumOut, pathOut, 
                     vDatumPath, subDomain)
        mesh = gpd.read_file(os.path.splitext(meshFile)[0]+'.shp', ignore_geometry = True) # not fully sure if it is the best way
    
    if epsgOut == 4326:
        sys.exit('Downscaling can not be done for the lat lon crs.')
    
    if clumpThreshold == 'from_mesh':
        thres = mesh.elemArea.min() * perMinElemArea
        unitEpsgOut = gdf.crs.axis_info[0].unit_name
        if 'metre' in unitEpsgOut or 'meter' in unitEpsgOut:
            pass
        else: ## feet
            thres = thres * 0.092903
    else:
        thres = clumpThreshold

    print(f'Static downscaling started')
    t1 = time.time()
    #print(f'    Shape file exported: {(t1 - t0)/60:0.2f} min')
    
    grassEnvVar(grassVer)
    ## import grass
    import grass.script as gs
    import grass.script.setup as gsetup
    
    if pathGrassLocation == None:
        pathGrassLocation = os.path.join(pathaux, 'grassLoc')
    
    print(f'    Start Setup grass environment')
    t11 = time.time()
    
    ## setup grass env
    setGrassEnv(grassVer, pathGrassLocation, createGrassLocation, gs, gsetup,
                pathRasFiles, rasterFiles, createLocMethod, epsgOut)
    
    t2 = time.time()
    print(f'    Setup grass environment: {(t2 - t11)/60:0.2f} min')
    
     ## setup growing
    print(f'    Start Downscaling preprocess')
    ## here the thres must be in square meters
    setupGrowing(pathOut, attrCol, exportMesh, meshFile, thres, gs)
    t3 = time.time()
    print(f'    Downscaling preprocess: {(t3 - t2)/60:0.3f} min')
    
    ## grow
    print(f'    Start growing')
    staticGrow(repLenGrowing, gs, os.path.splitext(os.path.basename(meshFile))[0])
    t4 = time.time()
    print(f'    Ready with static grow: {(t4 - t3)/60:0.3f} min')
    
    
   ## postprocess
    print(f'    Start postprocessing')
    postProcessStatic(compAdcirc2dem, floodDepth, pathOut, thres, gs, ras2vec)
    t5 = time.time()
    print(f'    Ready with postprocess: {(t5 - t4)/60:0.3f} min')
    print(f'Ready with static downscaling: {(t5 - t1)/60:0.3f} min')
    print(f'Kalpana finished sucsesfully after: {(t5 - t0)/60:0.3f} min')
    print(f'Output files saved on: {pathaux}')