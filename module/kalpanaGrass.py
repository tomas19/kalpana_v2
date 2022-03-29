import shutil
import os
import subprocess
import sys

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
    else:
        print("Warning: if DEM raster resolutions do not match, the aggregate DEM " \
                "resolution will match the resolution of the first input raster.")
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

def createGrassLoc(grassVer, locPath, createLocMethod='from_epsg', myepsg=4326, rasFile=None):
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
    existingLoc = locPath
    if os.path.isdir(existingLoc):
        print(f'Location: {existingLoc} will be deleted')
        shutil.rmtree(existingLoc)
    else:
        print('New location will be created')

    #### epsg can be obtained from a raster or specified by the user
    if createLocMethod == 'from_epsg':
        startCmd = grassBin + ' -c epsg:' + str(myepsg) + ' -e ' + locPath
    elif createLocMethod == 'from_raster':
        print(f'Projection from {rasFile}')
        startCmd = grassBin + ' -c ' + rasFile + ' -e ' + locPath
    else:
        sys.exit('The create location method specified is incorrect. See the docstring!')

    p = subprocess.Popen(startCmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = p.communicate()

    if p.returncode != 0:
        print(f'ERROR: {err}', file = sys.stderror)
        print(f'"ERROR: Cannot create location {startCmd}', file = sys.stderror)
        sys.exit(-1)
    else:
        print(f'Created location {locPath}')

def initGrass(locPath, mapset='PERMANENT'):
    ''' Set grass working environment
        Parameters
            locPath: str
                path of the grass location
            mapset: str
                mapset where the core data of the project can be stored
    '''
    #### import grass packages
    import grass.script.setup as gsetup
    #### init grass environment
    gsetup.init(f'{locPath}\\{mapset}')

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

def createGrassRegion(rasList, pkg):
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
        pkg.run_command('g.region',
                       raster=rasList)#,
    #                    quiet=True)
        ### patch rasters
        pkg.run_command('r.patch', 
                        input=rasList,
                        output='dem',
                        overwrite=True,
                        quiet=True)
    else:
        pkg.run_command('g.region', 
                       raster=rasList, 
                       quiet=True)
        #Set previous loop results to oldCostMap
        pkg.run_command('g.rename', 
                        raster=f'{rasList[0]},dem', 
                        overwrite = True,
                        quiet = True)
                        
                        

def vertUnitConvert(conv):
    ''' Convert vertical units from meter 2 feet or feet 2 meter.
        Parameters
            conv: str
                "m2ft" or "ft2m"
        Return
            None
    '''
    #import grass.script as gs
    gs.run_command('g.rename',raster="dem,demPreConv")
    if conv == 'm2ft':
        gs.mapcalc("$output=if(!isnull($demPreConv),$demPreConv*3.2808399,null())",
                        output="dem",
                        demPreConv="demPreConv@PERMANENT",
                        overwrite=True,
                        quiet=True)
    elif conv == 'ft2m':
        gs.mapcalc("$output=if(!isnull($demPreConv),$demPreConv/3.2808399,null())",
                        output="dem",
                        demPreConv="demPreConv@PERMANENT",
                        overwrite=True,
                        quiet=True)
    else:
        sys.exit('Wrong argument. Available options are: m2ft or ft2m')
    gs.run_command('g.remove',flags="f",type="all",name="demPreConv")
