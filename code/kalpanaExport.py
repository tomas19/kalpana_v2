import fiona
import xarray as xr
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import geopandas as gpd
import pandas as pd
from shapely.geometry import Polygon, LineString
from shapely import ops
from tqdm import tqdm
import itertools
import netCDF4 as netcdf
import simplekml

def contours2gpd(ncObj, var, levels, epsg, its=None, subDom=None):
    ''' Dataset to GeoDataFrame with contours as shapely LineStrings
        Parameters
            ncObj: xarray dataset
                adcirc input file
            var: string
                name of the variable to export
            levels: np.array
                contour levels
            epsg: int
                coordinate system
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
            indexSux: np.array
                array with indices of the nodes inside the requested subdomain.
                Output of the ncSubset function. If None (default), all domain is exported.
        Returns
            gdf: GeoDataFrame
                Linestrings as geometry and contour value in the "value" column.
                The name of the variable, the date and other info are also included in the columns
    '''
    nv = ncObj['element'][:,:] - 1 ## triangles starts from 1
    x = ncObj['x'][:].data
    y = ncObj['y'][:].data
    
    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    
    ## get contourlines
    if its == None: ## not time-varying file
        aux = ncObj[var][0, :].to_masked_array()

    else:
        aux = ncObj[var][its, :].to_masked_array()
        
    contours = plt.tricontour(tri, aux, levels = levels)
    plt.close()
    
    ## iteration over lines
    geoms, vals = [], []
    for icon, con in enumerate(contours.collections):
        val = contours.levels[icon]
        paths = con.get_paths()
        ## dismiss lines with less than 2 vertices
        aux0 = [LineString(path.vertices) for path in paths if len(path.vertices) > 1]
        aux1 = [val]*len(aux0)
        geoms.append(aux0)
        vals.append(aux1)
    
    ## list of lists to one list
    lines = list(itertools.chain(*geoms))
    values = list(itertools.chain(*vals))
    
    ## define GeoDataFrame
    gdf = gpd.GeoDataFrame(crs = f'EPSG:{epsg}', geometry = lines, index = range(len(lines)))
    gdf['value'] = values
    gdf['variable'] = [ncObj[var].name]*len(lines)
    units = ncObj[var].units
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in values]
    
    if its != None:
        t0 = ncObj['time'][0]
        ti = ncObj['time'][its]
        
        gdf['nTimeStep'] = [its]*len(lines)
        gdf['date'] = [str(ti.data).split('.')[0].replace('T', ' ')]*len(lines)
        gdf['nHours'] = [(ti.data - t0.data).item()/1e9]*len(lines)
        
    if subDom is not None:
        subDom = readSubDomain(subDom, epsg)
        gdf = gpd.clip(gdf, subDom)
    
    return gdf

def classifyPolygons(polys):
    ''' Classify polygons based on signed area to get inner and outer polygons.
        Used in filledContours2gpd.
        Parameters
            polys: list
                list with pair of coordinates of the polygon vertices
        Returns
            outer, inner: list
                list with vertices of each type of polygon
    '''
    outer = []
    inner = []
    areas = []
    for p in polys:
        area = signedArea_fx(p)
        areas.append(area)
        #print(area)
        if area >= 0:       
            outer.append(p)
        else:
            inner.append(p)
    
    return outer, inner

def signedArea_fx(ring):
    ''' Return the signed area enclosed by a ring in linear time using the algorith at
        https://web.archive.org/web/20080209143651/http://cgafaq.info:80/wiki/Polygon_Area.
        Used in classifyPolygons.  
        Parameters
            ring: list
                list with pairs of coordinates as tuples
        Returns
            signedArea: float
                area calculated using the vertices coordinates
    '''
    
    v2 = np.roll(ring, -1, axis=0)
    signedArea = np.cross(ring, v2).sum() / 2.0
    
    return signedArea

def pointsInsidePoly(points, polygon):
    ''' Get the subset of a cloud of points which are inside a polygon
        Used in filledContours2gpd.
        Parameters
            points: list
                list of pair of coordinates as tuples
            polygon: list
                list of coordinates of the polygon vertices as tuples
        Returns
            cont: array
                array with booleans
    '''
    p = mpl.path.Path(polygon)
    cont = p.contains_points(points)
    return cont

def filledContours2gpd(ncObj, var, levels, epsg, its=None, subDom=None):
    ''' Dataset to GeoDataFrame with filled contours as shapely polygons.
        TODO: WRITE THIS FX IN A MORE EASY-TO-READ WAY
        Parameters
            ncObj: xarray dataset
                Adcirc input file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            epsg: int
                coordinate system            
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
            indexSux: np.array
                array with indices of the nodes inside the requested subdomain.
                Output of the ncSubset function. If None (default), all domain is exported.
        Returns
            gdf: GeoDataFrame
                Polygons as geometry and contours min, average, max values as columns. 
                The name of the variable, the date and other info are also included in the columns
    '''
    
    nv = ncObj['element'][:,:] - 1 ## triangles starts from 1
    x = ncObj['x'][:].data
    y = ncObj['y'][:].data
    
    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    
    ## add max value over time and domain to the levels, just for polting porpuses.
    maxmax = np.max(ncObj[var])
    if maxmax > levels[-1]:
        levels = np.append(levels, [np.ceil(maxmax)])
    
    ## get contourlines
    if its == None: ## not time-varying file
        aux = ncObj[var][0, :].to_masked_array()

    else:
        aux = ncObj[var][its, :].to_masked_array()
    
    # aux2 = np.copy(aux)
    aux.data[np.isnan(aux.data)] = -99999.0
    
    contoursf = plt.tricontourf(tri, aux.data, levels = levels)
    plt.close()  
    
    geoms = []
    for icoll, coll in enumerate(contoursf.collections):
        ## get min an max malue of the interval
        vmin, vmax = contoursf.levels[icoll:icoll+2]
        vmean = np.mean([vmin, vmax])
        
        for p in coll.get_paths():
            p.simplify_threshold = 0.0
            # Removing polygons with less than 3 vertices
            polys = [g for g in p.to_polygons() if g.shape[0] >= 3] 
            # classify polygons based on the signed area
            outer,inner = classifyPolygons(polys)
            if len(inner) > 0:
                ## get first point of the inner polygon
                inner_points = [pts[0] for pts in inner]
            ## array of booleans
            overall_inout = np.zeros((len(inner),), dtype = bool)

            for out in outer:
                if len(inner) > 0:
                    ## check which inner polygons are inside out
                    inout = pointsInsidePoly(inner_points, out)
                    ## update boolean array
                    overall_inout = np.logical_or(overall_inout, inout)
                    ## use boolean array to discard inner polygons
                    out_inner = [g for f, g in enumerate(inner) if inout[f]]
                    ## remove holes or inner polygons
                    poly = Polygon(out, out_inner)
                else: ## no inner polygons
                    poly = Polygon(out)

                # clean-up polygons (remove intersections)
                if not poly.is_valid:
                    poly = poly.buffer(0.0)
                if poly.is_empty:
                    continue

                geoms.append((poly, vmin, vmax, vmean))
            
            # collect all interiors which do not belong to any of the exteriors
            outer_interiors = [interior for s, interior in enumerate(inner) if not overall_inout[s]]
            for k in outer_interiors:
                poly = Polygon(k[::-1]) ## reverse geometrty
                # clean-up polygons (remove intersections)
                if not poly.is_valid:
                    poly = poly.buffer(0.0)
                if poly.is_empty:
                    continue
                geoms.append((poly, vmin, vmax, vmean))
    
    ## define geopandas
    data = list(zip(*geoms))
    gdf = gpd.GeoDataFrame(crs = f'EPSG:{epsg}', geometry = list(data[0]), index = range(len(data[0])))
    # gdf['minVal'] = data[1]
    gdf['avgVal'] = data[3]
    # gdf['maxVal'] = data[2]
    gdf['variable'] = [ncObj[var].name]*len(data[3])
    units = ncObj[var].units
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in data[3]]
    
    if its != None:
        t0 = ncObj['time'][0]
        ti = ncObj['time'][its]
        
        gdf['nTimeStep'] = [its]*len(data[3])
        gdf['date'] = [str(ti.data).split('.')[0].replace('T', ' ')]*len(data[3])
        gdf['nHours'] = [(ti.data - t0.data).item()/1e9]*len(data[3])
        
    if subDom is not None:
        subDom = readSubDomain(subDom, epsg)
        gdf = gpd.clip(gdf, subDom)
    
    return gdf

def checkTimeVarying(ncObj):
    ''' Check if an adcirc input is time-varying or not.
        Parameters
            ncObj: xarray dataset
                Adcirc input file
        Returns
            timeVar: int
                1 if time-varying, 0 if not
    '''
    if ncObj['time'].shape[0] <= 1:
        ## not time-varying
        timeVar = 0
    elif (ncObj['time'][-1].data - ncObj['time'][0].data).astype(int) == 0:
        ## time variable has lenght 2 but dates are the same --> not time-varying
        timeVar = 0
    else:
        ## time-varying file
        timeVar = 1
    
    return timeVar

def runExtractContours(ncObj, var, levels, conType, epsg, subDom=None):
    ''' Run "contours2gpd" or "filledContours2gpd"
        Parameters
            ncObj: xarray dataset
                Adcirc input file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
            indexSux: np.array
                array with indices of the nodes inside the requested subdomain.
                Output of the ncSubset function. If None (default), all domain is exported.
        Returns
            gdf: GeoDataFrame
                Polygons or polylines as geometry columns. If the requested file is time-varying the GeoDataFrame will include all timesteps.
            
    '''
    timeVar = checkTimeVarying(ncObj)
    if timeVar == 0:
        if conType == 'polyline':
            gdf = contours2gpd(ncObj, var, levels, epsg, subDom = subDom)
        elif conType == 'polygon':
            gdf = filledContours2gpd(ncObj, var, levels, epsg, subDom = subDom)
        else:
            print('only "polyline" and "polygon" types are supported!')
            sys.exit(-1)
    else:
        listGdf = []
        if conType == 'polyline':
            for t in tqdm(range(ncObj['time'].shape[0])):
                aux = contours2gpd(ncObj, var, levels, epsg, t, subDom)
                listGdf.append(aux)
        elif conType == 'polygon':
            for t in tqdm(range(ncObj['time'].shape[0])):
                aux = filledContours2gpd(ncObj, var, levels, epsg, t, subDom)
                listGdf.append(aux)
        else:
            print('only "polyline" and "polygon" types are supported!')
            sys.exit(-1)
        gdf = gpd.GeoDataFrame(pd.concat(listGdf, axis=0, ignore_index=True), crs=epsg)
    
    return gdf

def mesh2gdf(ncObj, epsg):
    ''' Write adcirc mesh as GeoDataFrame and extract centroid of each element. Used to create submesh
        Parameters:
            ncObj: xarray dataset
                adcirc file
            epsg: int
                coordinate system
        Returns
            gdf: GeoDataFrame
                GeoDataFrame with polygons as geometry and centroid coordinates in columns
        
    '''
    nv = ncObj['element'][:,:] - 1 ## triangles starts from 1
    x = ncObj['x'][:].data
    y = ncObj['y'][:].data
    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    xvertices = x[tri.triangles[:]]
    yvertices = y[tri.triangles[:]]
    listElem = np.stack((xvertices, yvertices), axis = 2)
    pols = [Polygon(x) for x in listElem]
    gdf = gpd.GeoDataFrame(geometry = pols, crs = f'EPSG:{epsg}')
    gdf['centX'] = xvertices.mean(axis = 1)
    gdf['centY'] = yvertices.mean(axis = 1)
    
    return gdf
    
def readSubDomain(subDomain, epsg):
    ''' Read a user specified subdomain and transform it to GeoDataFrame
        Parameters
            subDomain: str or list
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            epsg: int
                epsg code of the used system of coordinates       
        Returns
            gdfSubDomain: GeoDataFrame
                user specified subdomain as geodataframe
    '''
    if type(subDomain) == str:
        ## kml or shp file
        if subDomain.endswith('.shp'):
            gdfSubDomain = gpd.read_file(subDomain)
        elif subDomain.endswith('.kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            gdfSubDomain = gpd.read_file(subDomain, driver = 'KML')
        else:
            print('Only shape and kml format are suported for sub domain generation!')
            sys.exit(-1)
        xAux, yAux = gdfSubDomain.geometry[0].exterior.coords.xy
        extCoords = list(zip(xAux, yAux))
    
    elif type(subDomain) == list and len(subDomain) == 4:
        ## only UL lon, UL lat, LR lon LR lat
        ulLon, ulLat, lrLon, lrLat = subDomain
        extCoords = [(ulLon, ulLat), (ulLon, lrLat), (lrLon, lrLat), (lrLon, ulLat), (ulLon, ulLat)]
        poly = Polygon(extCoords)
        gdfSubDomain = gpd.GeoDataFrame(geometry = [poly], crs = f'EPSG:{epsg}')
    else:
        print('subDomain must be the path of a kml or shapefile, or a list with the coordinates of ' \
              'the upper left and lower right corners of a box')
        sys.exit(-1)
    
    return gdfSubDomain
    
def ncSubset(ncObj, subDomain, epsg):
    ''' Generate a subdomain of the original adcirc input mesh based on a user input.
        Parameters
            ncObj: xarray dataset
                adcirc input file
            subDomain: str or list
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            epsg: int
                epsg code of the used system of coordinates
        Returns
            nodesInside: array
                indices of the nodes within the requested subdomain.
            elemInside: array
                indices of the elements within the requested subdomain.
    '''
    subDomain = readSubDomain(subDomain, epsg)
    xAux, yAux = subDomain.geometry[0].exterior.coords.xy
    extCoords = list(zip(xAux, yAux))
    
    x = ncObj['x'][:].data
    y = ncObj['y'][:].data
    nodes = list(zip(x, y))
    
    inside = pointsInsidePoly(nodes, extCoords)
    nodesInside = np.where(inside)[0]
    
    mesh = mesh2gdf(ncObj, epsg)    
    elemCent = list(zip(mesh.centX, mesh.centY))
    boolInside = pointsInsidePoly(elemCent, extCoords)
    elemInside = np.where(boolInside)[0]    
    
    return nodesInside, elemInside
    
def meshSubset(ncObj, subDomain, epsg):
    ''' Generate a subdomain of the original adcirc input mesh based on a subdomain
        Parameters
            mesh: GeoDataFrame
                output of the mesh2gdf function
            subDomain: str or list
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            epsg: int
                epsg code of the used system of coordinates
        Returns
            gdfMeshSubset: GeoDataFrame
                GeoDataFrame with the mesh subset
    '''
    
    mesh = mesh2gdf(ncObj, epsg)
    
    sub = readSubDomain(subDomain, epsg)
    xAux, yAux = sub.geometry[0].exterior.coords.xy
    extCoords = list(zip(xAux, yAux))
    
    elemCent = list(zip(mesh.centX, mesh.centY))
    boolInside = pointsInsidePoly(elemCent, extCoords)
    elemInsideIx = np.where(boolInside)[0]
    gdfMeshSubset = mesh.iloc[elemInsideIx, :]
    
    return gdfMeshSubset
    
def nc2xr(ncFile, var):
    ''' Write netcdf as xarray dataset. May be redundant but I had some problems reading adcirc files with xarray.
        Xarray is faster and easier to work with than netCDF4 object.
        Parameters
            ncFile: string
                complete path of the adcirc input
            var: string
                name of the variable to work with
        Returns
            ds: xarray dataset
                netCDF4 object casted as xarray dataset
    '''
    with netcdf.Dataset(ncFile, 'r') as ncObj:
    
        epoch = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        dates = [epoch + pd.Timedelta(seconds = float(x)) for x in ncObj['time'][:]]
    
        if checkTimeVarying(ncObj) == 1:
            dummy1 = ncObj[var][:,:].data
            nnodes = range(ncObj[var].shape[1])
            dims = ['time', 'node']
            # coords = {'time': ncObj['time'][:].data, 'node': nnodes}
            coords = {'time': dates, 'node': nnodes}
        else:
            dummy1 = ncObj[var][:].data
            dummy1 = dummy1.reshape((1, len(dummy1)))
            nnodes = range(ncObj[var].shape[0])
            dims = ['node']
            coords = {'time': [0], 'node': nnodes}            
        
        dummy2 = np.ma.masked_invalid(dummy1)
        dummy3 = dummy2.filled(fill_value = -99999.0)
        
        ds = xr.Dataset({
                        var:xr.DataArray(data = dummy3,
                                          dims = ['time', 'node'],
                                          coords = coords),
                       'element': xr.DataArray(data = ncObj['element'][:, :].data,
                                              dims = ['nele', 'nvertex']),
                       'x': xr.DataArray(data = ncObj['x'][:].data,
                                         dims = ['node'],
                                         coords = {'node': nnodes}),
                       'y': xr.DataArray(data = ncObj['y'][:].data,
                                         dims = ['node'],
                                         coords = {'node': nnodes})                                     
                        })
        
        ds[var].attrs['name'] = ncObj[var].long_name.capitalize()
        ds[var].attrs['units'] = ncObj[var].units
        # ds['time'].attrs['base_date'] = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        # ds['time'].attrs['units'] = f'Seconds since {epoch}'
        
    return ds
    
def nc2shp(ncFile, var, levels, conType, epsg, pathOut, subDomain=None, dateSubset=None):
    '''
    '''
    nc = nc2xr(ncFile, var)
    levels = np.arange(levels[0], levels[1], levels[2])
    if subDomain == None:
        gdf = runExtractContours(nc, var, levels, conType, epsg)
    else:
        # indexSub = ncSubset(nc, var, subDomain, epsg)
        gdf = runExtractContours(nc, var, levels, conType, epsg, subDomain)
    
    
    if pathOut.endswith('.shp'):
        gdf.to_file(pathOut)
    elif pathOut.endswith('.gpkg'):
        gdf.to_file(pathOut, driver = 'GPKG')
    
    return gdf
    
####################################### GOOGLE EARTH ################################################

def createColorbar(levels, varName, units, cmap='viridis', fileName='tempColorbar.jpg', filePath='.'):
    ''' Create colorbar image to be imported in the kml file
        Parameters
            levels: np.array
                contour levels
            varName: string
                name of the variable to export, just used for the legend.
            units: string
                units of the variable to export, just used for the legend.
            cmap: string
                colormap name
        Return
            None. The image will be dumped in the same path where the script is executed.
                
    '''
    norm = mpl.colors.Normalize(vmin = levels[0], vmax = levels[-1])
    cmap = mpl.cm.get_cmap(cmap)
    fig, ax = plt.subplots(figsize=(1, 8))
    cb = mpl.colorbar.ColorbarBase(ax, cmap = cmap, norm = norm, extend = 'neither',
                                    extendfrac = 'auto', ticks = levels, spacing = 'uniform',
                                    orientation = 'vertical')
    cb.set_label(f'{varName} [{units}]')
    fig.savefig(os.path.join(filePath, fileName), dpi = 300, bbox_inches = 'tight')
    plt.close()
    
def kmlScreenOverlays(kml, colorbar=True, colorbarFile='tempColorbar.jpg', logo=True, 
                        logoFile='logo.png', logoUnits='fraction', logoDims=None):
    ''' Create screen overlay for the kml file with the colorbar and the project logo.
        TODO: check logo chile path when calling the code from github repo
        Parameters
            kml: simplekml object.
                output of lines2kml or poly2kml functions.
            colorbar: bolean, default True.
                True for including the colorbar.
            logo: bolean, default True.
                True for including the logo.
            logoFile: string, default "logo.png"
                name of the file to be included as logo
            logoUnits: string, default fraction
                the other available option is pixel.
            logoDims: list, default None
                x and y dimensions of the logo overlay
        Returns
            None    
    '''
    if colorbar == True:
        screen1 = kml.newscreenoverlay(name='Colorbar')
        screen1.icon.href = 'tempColorbar.jpg'
        screen1.overlayxy = simplekml.OverlayXY(x= 0 , y = 0, xunits = simplekml.Units.fraction,
                                         yunits = simplekml.Units.fraction)
        screen1.screenxy = simplekml.ScreenXY(x = 0, y = 0.1, xunits = simplekml.Units.fraction,
                                         yunits = simplekml.Units.fraction)
        screen1.size.x = 0.08
        screen1.size.y = 0.65
        screen1.size.xunits = simplekml.Units.fraction
        screen1.size.yunits = simplekml.Units.fraction
    
    if logo == True:
        screen2 = kml.newscreenoverlay(name = 'logo')
        screen2.icon.href = logoFile
        screen2.overlayxy = simplekml.OverlayXY(x = 0, y = 1, xunits = simplekml.Units.fraction,
                                       yunits = simplekml.Units.fraction)
        screen2.screenxy = simplekml.ScreenXY(x = 0, y = 1, xunits = simplekml.Units.fraction,
                                        yunits = simplekml.Units.fraction)
        if logoUnits == 'fraction':     
            screen2.size.xunits = simplekml.Units.fraction
            screen2.size.yunits = simplekml.Units.fraction
        elif logoUnits == 'pixel':
            screen2.size.xunits = simplekml.Units.pixel
            screen2.size.yunits = simplekml.Units.pixel
        
        if logoDims == None:
            screen2.size.x = 0.85
            screen2.size.y = 0.08
        else:
            screen2.size.x = logoDims[0]
            screen2.size.y = logoDims[1]
            
def lines2kml(gdf, levels, cmap='viridis'):
    ''' Write GeoDataFrame with shapely linestrings as kml
        TODO: ADD SYS.EXIT IF GDF HAS A DATE COLUMN?
        Parameters
            gdf: GeoDataFrame
                output of runExtractContours function
            levels: np.array
                contour levels
            cmap: str, default "viridis"
                colormap
        Returns
            kml: simplekml object
    '''
    norm = mpl.colors.Normalize(vmin = levels[0], vmax = levels[-1])
    cmap = mpl.cm.get_cmap(cmap)
    m = mpl.cm.ScalarMappable(norm = norm, cmap = cmap)

    kml = simplekml.Kml()
    
    for i in tqdm(gdf.index):
        ls = kml.newlinestring(name = gdf.loc[i, 'labelCol'])
        try:
            coords = list(gdf.loc[i, 'geometry'].coords)
        except:
            aux = []
            for line in gdf.loc[i, 'geometry'].geoms:
                aux.extend(list(line.coords))
            # dummy = ops.linemerge(gdf.loc[i, 'geometry'])
            coords = aux
        ls.coords = coords
        value = gdf.loc[i, 'value']
        r, g, b, a = m.to_rgba(value)
        ls.style.linestyle.color = simplekml.Color.rgb(int(255*r), int(255*g), int(255*b))
        ls.style.linestyle.width = 2
        ls.description = gdf.loc[i, 'labelCol']
    
    return kml
    
def polys2kml(gdf, levels, cmap='viridis'):
    ''' Write GeoDataFrame with shapely polygons as kml
        Parameters
            gdf: GeoDataFrame
                output of runExtractContours function
            levels: np.array
                contour levels
            cmap: str, default "viridis"
                colormap
        Returns
            kml: simplekml object
    '''
    norm = mpl.colors.Normalize(vmin = levels[0], vmax = levels[-1])
    cmap = mpl.cm.get_cmap('viridis')
    m = mpl.cm.ScalarMappable(norm = norm, cmap = cmap)

    kml = simplekml.Kml()

    for i in tqdm(gdf.index):
        try:
            outerCoords = list(zip(gdf.loc[i, 'geometry'].exterior.coords.xy[0], 
                                        gdf.loc[i, 'geometry'].exterior.coords.xy[1]))

            innerCoords = list(gdf.loc[i, 'geometry'].interiors)    
            if len(innerCoords) > 0:
                innerCoords = [list(interior.coords) for interior in innerCoords]
            
            pol = kml.newpolygon(name = gdf.loc[i, 'labelCol'])
            pol.outerboundaryis = outerCoords
            pol.innerboundaryis = innerCoords
            value = gdf.loc[i, 'avgVal']
            r, g, b, a = m.to_rgba(value)
            col = simplekml.Color.rgb(int(255*r), int(255*g), int(255*b))
            pol.style.linestyle.color = col
            pol.style.linestyle.width = 2
            pol.description = gdf.loc[i, 'labelCol']
    #         pol.style.polystyle.color = col
            pol.style.polystyle.color = simplekml.Color.changealphaint(100, col)
        
        except:
            for poly in gdf.loc[i, 'geometry'].geoms:
                outerCoords = list(zip(poly.exterior.coords.xy[0], 
                                        poly.exterior.coords.xy[1]))

                innerCoords = list(poly.interiors)    
                if len(innerCoords) > 0:
                    innerCoords = [list(interior.coords) for interior in innerCoords]
                
                pol = kml.newpolygon(name = gdf.loc[i, 'labelCol'])
                pol.outerboundaryis = outerCoords
                pol.innerboundaryis = innerCoords
                value = gdf.loc[i, 'avgVal']
                r, g, b, a = m.to_rgba(value)
                col = simplekml.Color.rgb(int(255*r), int(255*g), int(255*b))
                pol.style.linestyle.color = col
                pol.style.linestyle.width = 2
                pol.description = gdf.loc[i, 'labelCol']
        #         pol.style.polystyle.color = col
                pol.style.polystyle.color = simplekml.Color.changealphaint(100, col)


    return kml
    
def nc2kmz(ncFile, var, levels, conType, epsg, pathOut, overlay=True, logoPath=None, subDomain=None, dateSubset=None, cmap='viridis'):
    '''
    '''
    nc = nc2xr(ncFile, var)
    
    if checkTimeVarying(nc) == 1:
        print('Time-varying files can not be exported as kmz!')
        sys.exit(-1)
    else:
        levels = np.arange(levels[0], levels[1], levels[2])
        
        if subDomain == None:
            gdf = runExtractContours(nc, var, levels, conType, epsg)
        else:
            # indexSub = ncSubset(nc, var, subDomain, epsg)
            gdf = runExtractContours(nc, var, levels, conType, epsg, subDomain)   
        
        if conType == 'polygon':
            kml = polys2kml(gdf, levels, cmap)
        elif conType == 'polyline':
            kml = lines2kml(gdf, levels, cmap)
        else:
            print('Only "polygon" o "polyline" formats are supported')
            sys.exit(-1)
        
        if overlay == True:
            name = nc[var].name
            units = nc[var].units
            createColorbar(levels, name, units)
            kmlScreenOverlays(kml, logoFile = logoPath)
            
        kml.savekmz(pathOut, format = False)
        os.remove('tempColorbar.jpg')

