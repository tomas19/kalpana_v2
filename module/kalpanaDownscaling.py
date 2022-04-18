import shutil
import os
import subprocess
import sys
import time
from kalpanaExport import nc2shp

'''
    All functions using GRASS GRIS has an argument 'pkg' which is the grass package.
    The package can't be imported before runing the 'grassEnvVar' function so it's
    not directly imported in this module and will be imported when the functions
    are executed.
'''

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
                epsg code (https://epsg.io/), only used if createLocMehot == "from_epsg". Default 4326.
            rasFiles: str
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
            mapset: str
                mapset where the core data of the project can be stored
    '''
    # #### import grass packages
    # import grass.script.setup as gsetup
    #### init grass environment
    pkg.init(f'{locPath}\\{mapset}')

def importRasters(rasFiles, rasterRes, pkg):
    ''' Import rasters using Grass gis.
        Parameters
            rasFiles: list
                list with complete path of each raster file
            rasterRes: str or int
                input "align" if to use the res of the FIRST raster in rasFiles list or
                an integer to use a user-specified value. TODO CHECK UNITS.
        Returns
            rasterOutList: list
                list with imported raster files
    '''
    #import grass.script as gs
    rasterOutList = []
    #### using first raster resolution
    if rasterRes == 'align':
        for ras in rasFiles:
            rasObj = os.path.basename(ras).split('.')[0]
            ### import
            pkg.run_command('r.import',
                            overwrite=True,
                            input=ras,
                            output=rasObj)
            rasterOutList.append(rasObj)

    elif type(rasterRes) == int:
        for ras in rasFiles:
            rasObj = os.path.basename(ras).split('.')[0]
            pkg.run_command('r.import',
                            overwrite=True,
                            input=ras,
                            output=rasObj,
                            resolution='value',
                            resolution_value=rasterRes,
                            resample='bilinear',
                            extent='input')
            rasterOutList.append(rasObj)

    else:
        print('Raster resolution options are "align" or an integer value')
        sys.exit(-1)
        
    return rasterOutList

def setDownscalingDEM(rasList, pkg):
    ''' Create grass region based in input rasters
        Parameters
            rasList: list
                list of the raster objects imported
        Return
            None
    '''
    #import grass.script as gs
    if len(rasList) > 1:
        ### set region
        # pkg.run_command('g.region',
                       # raster=rasList)#,
    #                    quiet=True)
        ### patch rasters
        pkg.run_command('r.patch', 
                        input=rasList,
                        output='dem',
                        overwrite=True,
                        quiet=True)
    else:
        # pkg.run_command('g.region', 
                       # raster=rasList, 
                       # quiet=True)
        #Set previous loop results to oldCostMap
        pkg.run_command('g.rename', 
                        raster=f'{rasList[0]},dem', 
                        overwrite = True,
                        quiet = True)
                        

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
                pathRasFiles, rasterFiles, createLocMethod, myepsg, rasterRes):
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
            rasterRes: str or int
                input "align" if to use the res of the FIRST raster in rasFiles list or
                an integer to use a user-specified value. TODO CHECK UNITS.
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
        grassRasList = importRasters(rasFiles, rasterRes, pkg0) ## try rasters with different resolution? check if rasters are in different crs that location
        print(f'        import raster: {(time.time() - ta)/60: 0.3f} min')
        
        ta = time.time()
        setDownscalingDEM(grassRasList, pkg0)
        print(f'        set downscaling dem: {(time.time() - ta)/60: 0.3f} min')
    
    else:
        ta = time.time()
        initGrass(pathGrassLocation, pkg1)
        print(f'        init grass: {(time.time() - ta)/60: 0.3f} min')
        
def setupGrowing(kalpanaShp, attrCol, pkg):
    ''' Preprocess kalpana shape file
        Parameters
            kalpanaShp: str
                path of the shape file with AdCirc results
            attrCol: str
                name of the attribute column
    '''
    ## import shape file with max water level
    # pkg.run_command('v.import', input = kalpanaShp, overwrite = True, 
                   # quiet = True, snap = 0.000001)
    pkg.run_command('v.in.ogr', input = kalpanaShp, overwrite = True,
                    quiet = True, snap = 0.000001, min_area = 10,
                    flags = 'o')
                    
    # reg = pkg.read_command('g.region', flags = 'g').split()
    
    # ewRes = float(reg[7].split("=")[1])
    # nsRes = float(reg[6].split("=")[1])
    # pkg.run_command('g.region', raster = 'dem', nsres = nsRes, 
                   # ewres = ewRes, overwrite = True, quiet = True)
                   
    pkg.run_command('v.to.rast', input = os.path.basename(kalpanaShp[:-4]), 
                    type = 'area', output = 'kalpanaRast', use = 'attr', 
                    quiet = True, attribute_column = attrCol, overwrite = True)
    # pkg.run_command('r.mask', raster = 'dem', quiet = True, 
                   # overwrite = True)
                   
def staticGrow(growRadius, pkg):
    ''' Execture r.grow.distance algorithm
        Parameters
            growRadius: int or float
                growing radius must be specified in meters. If 0 the grown distance is not limited to the
                growRadius, so results will be downscaled until the land level is higher than the AdCirc removing
                isolated cells. If positive, results will be downscaled not far way of the growRadius. If negative,
                results will be shrink by growRadius.
    '''
    if growRadius == 0: ## grow raster cells without a limit distance
        pkg.run_command('r.grow.distance', input = 'kalpanaRast', value = 'grownRastVal', overwrite = True,
                        quiet = True)
        
        pkg.mapcalc("$output = if(!isnull($input), $input, if($base < $new, $new, null()))",
                    output = 'grownKalpanaRast0', input = 'kalpanaRast', base = 'dem', 
                    new = 'grownRastVal', quiet = True, overwrite = True)
        
    elif growRadius > 0: ## grow raster cells with a limiting distance on growRadius
        growRadiusSq = growRadius**2 #distance on r.grow is calculated using squared metric so it is the squared of the actual distance
        pkg.run_command('r.grow.distance', input = 'kalpanaRast', metric = 'squared', distance = 'grownRastDist',
                        value = 'grownRastVal', overwrite = True, flags = 'm', quiet = True) ## flag m: distance in meters
                      
        pkg.mapcalc("$output = if(!isnull($input), $input, if($dist < $radius && $base < $new, $new, null()))",
                  output = 'grownKalpanaRast0', input = 'kalpanaRast', radius = growRadiusSq, base = 'dem',
                  new = 'grownRastVal', dist = 'grownRastDist', quiet = True, overwrite = True)
    
    else: #growradius < 0 --> shrink results, not sure why it is usefull but it was in the original version of kalpana
        kv = pkg.region()
        scale = np.sqrt(kv['nsres'] * kv['ewres'])
        growRadiusPixelsSq = (growRadius * scale)**2
        pkg.run_command('r.grow.distance', input = 'kalpanaRast', metric = 'squared', distance = 'grownRastDist', 
                          value = 'grownRastVal', flags = 'n', quiet = True,
                          overwrite = True)  ## flag n: distance to the nearest NO-NULL cell. value in pixels
        
        pkg.mapcalc("$output = if($dist < $radius, null(), $input)", output = 'grownKalpanaRast0', 
                    radius = growRadiusPixelsSq, input = 'kalpanaRast', dist = 'grownRastDist', 
                    quiet = True, overwrite = True)
                   
def postProcessStatic(compAdcirc2dem, floodDepth, kalpanaShp, grownRadius, pkg, clumpSizeThreshold = 50):
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
            growRadius: int or float
                growing radius must be specified in meters. It is used to save the
                downscaled shp and raster files
            threshold: int. Default 50
                threshold to delete isolated cells. Groups of cells smaller than threshold
                will be removed if are not connected to the main water area.
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
    # Gives all non-null cells in grown raster a single uniform value (-1).
    pkg.mapcalc("$output = if(!isnull($input), -1, null())", output = 'temp1', 
               input = 'grownKalpanaRast1', quiet = True, overwrite = True)
    print(f'        Assign -1: {(time.time() - ta)/60:0.3f}')
    
    #Groups the uniform raster by giving each connected group of cells a unique ID.
    #The goal is to remove isolated clumps not connected to original raster.
    ta = time.time()
    pkg.run_command('r.clump', input = 'temp1', output = 'temp2', 
                   quiet = True, overwrite = True)
    print(f'        Clump: {(time.time() - ta)/60:0.3f}')
    
    # Identifies original clumps found in the ADCIRC raster.
    ta = time.time()
    pkg.mapcalc("$output = if(!isnull($A) && !isnull($B), $B, null())", output = 'temp3', 
               A = 'kalpanaRast', B = 'temp2', quiet = True, overwrite = True)
    print(f'        Identifies original clumps found in the ADCIRC raster: {(time.time() - ta)/60:0.3f}')
    
    # Sorts clump areas largest to smallest, removes largest clump (null clump).
    # List format: [clump1, #cells1, clump2, #cells2, ...].
    ta = time.time()
    areas = pkg.read_command('r.stats', input = 'temp3', sort = 'desc', 
                            flags = 'c', quiet = True).split()
    
    ncells = [x for x in areas[3::2] if int(x) >= 50]
    clumpsID = [x for x, y in zip(areas[2::2], areas[1::2]) if int(y) >=  clumpSizeThreshold]
    #Interleave the two lists 
    areas = [val for tup in zip(clumpsID, ncells) for val in tup]
    print(f'        find areas: {(time.time() - ta)/60:0.3f}')
    
    # Assigns uniform value (-1) to all areas that include and are connected to
    # original ADCIRC raster (r.reclass requires a list of rules as str)
    #Format the large area IDs for use with r.reclass
    ta = time.time()
    reclassList = ''.join([f"{i} = -1\n" for i in areas[::2]])
    pkg.write_command('r.reclass', input = 'temp2', output = 'temp4', rules = '-', 
                stdin = reclassList, quiet = True, overwrite = True)
    print(f'        Reclass: {(time.time() - ta)/60:0.3f}')
    
    ta = time.time()
    # Passes back grown ADCIRC cells if they coincide with the assigned value (-1).
    pkg.mapcalc("$output = if($A == -1, $B, null())", output = 'grownKalpanaRastLevel', A = 'temp4', 
               B = 'grownKalpanaRast1', quiet = True, overwrite = True)
    print(f'        Passes back grown ADCIRC: {(time.time() - ta)/60:0.3f}')
    
    
    pathOut = os.path.dirname(kalpanaShp)
    fileOut = os.path.basename(kalpanaShp)[:-4]# remove extension
    
    ta = time.time()
    # export raster as tif
    pkg.run_command('r.out.gdal', input = 'grownKalpanaRastLevel', flags = 'm', format = 'GTiff', nodata = -9999, 
                   output = os.path.join(pathOut, f'{fileOut}_level_downscaled_Radius{grownRadius}m.tif'), 
                   overwrite = True)
    print(f'        export as tif level: {(time.time() - ta)/60:0.3f}')
    
    ta = time.time()
    #Export to ESRI shapefile
    pkg.run_command('r.to.vect', input = 'grownKalpanaRastLevel', output = 'grownKalpanaVectLevel', type = 'area')
    print(f'        ras to vector: {(time.time() - ta)/60:0.3f}')
    
    ta = time.time()
    pkg.run_command('v.out.ogr', input = 'grownKalpanaVectLevel', 
                      output = os.path.join(pathOut, f'{fileOut}_level_downscaled_Radius{grownRadius}m'), 
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
                       output = os.path.join(pathOut, f'{fileOut}_depth_downscaled_Radius{grownRadius}m.tif'), 
                       overwrite = True)
        print(f'        export as tif depth: {(time.time() - ta)/60:0.3f}')
        
        ta = time.time()
        #Export to ESRI shapefile
        pkg.run_command('r.to.vect', input = 'grownKalpanaRastDepth', output = 'grownKalpanaVectDepth', type = 'area')
        print(f'        ras to shp depth: {(time.time() - ta)/60:0.3f}')
        
        ta = time.time()
        pkg.run_command('v.out.ogr', input = 'grownKalpanaVectDepth', 
                        output = os.path.join(pathOut, f'{fileOut}_depth_downscaled_Radius{grownRadius}m'), 
                        type = 'area', format = 'ESRI_Shapefile', flags = 'se', quiet = True, overwrite = True)
        print(f'        export as shp depth: {(time.time() - ta)/60:0.3f}')
        
def runStatic(ncFile, levels, epsgIn, epsgOut, pathOut,  grassVer, pathRasFiles, rasterFiles,
              var = 'zeta_max', conType = 'polygon', npro = 1, subDomain = None, pathGrassLocation = None, createGrassLocation = True,
              createLocMethod = 'from_raster', rasterRes = 'align', attrCol = 'avgVal', growRadius = 500, compAdcirc2dem = True, 
              floodDepth = True, clumpThreshold = 50):
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
            pkg0: grass gis module. DEFAULT gs
                gs from running import grass.script as gs
            pkg1: grass gis module. DEFAULT gsetup
                gsetup from running import grass.script.setup as gsetup
        ********************************************************************************************************************
        ***************************************** OPTIONAL inputs of nc2shp function ***************************************
        ********************************************************************************************************************
            var: string. DEFAULT zeta_max
                Name of the variable to export
            conType: string. DEFAULT polygon
                'polyline' or 'polygon'
            npro: int. DEFAULT 1 (running in serial)
                number of worker processes. More info: https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.
                For using all available processors, input 999.
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates. The crs must be the same of the
                adcirc input file.
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
            rasterRes: str or int
                input "align" if to use the res of the FIRST raster in rasFiles list or
                an integer to use a user-specified value. TODO CHECK UNITS.
            attrCol: str. DEFAULT 'avgVal'
                name of the attribute column
            growRadius: int or float. DEFAULT 50
                growing radius must be specified in meters. If 0 the grown distance is not limited to the
                growRadius, so results will be downscaled until the land level is higher than the AdCirc levels. 
                If positive, results will be downscaled limitinf the distance by growRadius. If negative,
                results will be shrink by growRadius.
            compAdcirc2dem: boolean. DEFAULT True
                True for removing ADCIRC cells with values lower than the dem.
            flooDepth: boolean . DEFAULT True
                True for transform water levels to water depth. False for export water levelss.
            clumpThreshold: int. Default 50
                threshold to delete isolated cells. Groups of cells smaller than threshold
                will be removed if are not connected to the main water area.
    '''
    t0 = time.time()
    print(f'Static downscaling started')
    pathaux = os.path.dirname(pathOut)
    if not os.path.exists(pathaux):
        os.mkdir(pathaux)

    gdf = nc2shp(ncFile, var, levels, conType, epsgIn, epsgOut, pathOut, npro, subDomain)
    t1 = time.time()
    print(f'    Shape file exported: {(t1 - t0)/60:0.2f} min')
    
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
                pathRasFiles, rasterFiles, createLocMethod, epsgOut, rasterRes)
    t2 = time.time()
    print(f'    Setup grass environment: {(t2 - t11)/60:0.2f} min')
    
     ## setup growing
    print(f'    Start Downscaling preprocess')
    setupGrowing(pathOut, attrCol, gs)
    t3 = time.time()
    print(f'    Downscaling preprocess: {(t3 - t2)/60:0.3f} min')
    
    ## grow
    print(f'    Start growing')
    staticGrow(growRadius, gs)
    t4 = time.time()
    print(f'    Ready with static grow: {(t4 - t3)/60:0.3f} min')
    
   ## postprocess
    print(f'    Start postprocessing')
    postProcessStatic(compAdcirc2dem, floodDepth, pathOut, growRadius, gs, clumpThreshold)
    t5 = time.time()
    print(f'    Ready with postprocess: {(t5 - t4)/60:0.3f} min')
    print(f'Ready with static downscaling: {(t5 - t0)/60:0.3f} min')
    print(f'Output files saved on: {pathaux}')