import fiona
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
from multiprocessing import Pool

def contours2gpd(ncObj, var, levels, epsg, its, subDom):
    ''' Dataset to GeoDataFrame with contours as shapely LineStrings
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
                adcirc input file
            var: string
                name of the variable to export
            levels: np.array
                contour levels
            epsg: int
                coordinate system
            its: int
                timestep to be exported. If None it is assumed that the file is not time-varying.
            subDom: str or list.
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
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
        aux = ncObj[var][:]#.to_masked_array()

    else:
        aux = ncObj[var][its, :]#.to_masked_array()
        
    contours = plt.tricontour(tri, aux.data, levels = levels)
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
    gdf['name'] = [ncObj[var].long_name.capitalize()]*len(lines)
    units = ncObj[var].units
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in values]
    
    if its != None:
        # t0 = ncObj['time'][0]
        t0 = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        # ti = ncObj['time'][its]
        ti = pd.Timedelta(seconds = int(ncObj['time'][its]))
        
        gdf['nTimeStep'] = [its]*len(lines)
        # gdf['date'] = [str(ti.data).split('.')[0].replace('T', ' ')]*len(lines)
        gdf['date'] = [str(t0 + ti)]*len(lines)
        # gdf['nHours'] = [(ti.data - t0.data).item()/1e9]*len(lines)
        gdf['nHours'] = [ti.total_seconds()/3600]*len(lines)
        
    if subDom is not None:
        subDom = readSubDomain(subDom, epsg)
        gdf = gpd.clip(gdf, subDom)
    
    return gdf
    
def contours2gpd_mp(data, x, y, nv, name, units, longName, epoch, times, levels, epsg, subDom, its):
    ''' Run contours2gpd function with multiprocessing. NetCDF object can not be used with the Pool class,
        that's the reasson why only arrays are used as arguments.
        Parameters
            data: array
                values of the adcirc input file
            x: array
                longitude coordinates
            y: array
                latitude coordinates
            nv: array
                mesh elements
            name: string
                variable long name
            epoch: timestamp
                start date
            times: array
                time vector
            levels: np.array
                contour levels
            epsg: int
                coordinate system
            subDom: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
        Returns
            gdf: GeoDataFrame
                Linestrings as geometry and contour value in the "value" column.
                The name of the variable, the date and other info are also included in the columns
    '''
    nv = nv - 1 ## triangles starts from 1
    
    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    
    ## get contourlines
    aux = data[its, :]
        
    contours = plt.tricontour(tri, aux.data, levels = levels)
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
    gdf['variable'] = [name]*len(lines)
    gdf['name'] = [longName.capitalize()]*len(lines)
    units = units
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in values]
    
    ti = pd.Timedelta(seconds = times[its])
        
    gdf['nTimeStep'] = [its]*len(lines)
    gdf['date'] = [str(epoch + ti)]*len(lines)
    gdf['nHours'] = [ti.total_seconds()/3600]*len(lines)
        
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

def filledContours2gpd(ncObj, var, levels, epsg, its, subDom):
    ''' Dataset to GeoDataFrame with filled contours as shapely polygons.
        TODO: WRITE THIS FX IN A MORE EASY-TO-READ WAY
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
                Adcirc input file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            epsg: int
                coordinate system            
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
            subDom: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
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
        aux = ncObj[var][:]#.to_masked_array()

    else:
        aux = ncObj[var][its, :]#.to_masked_array()
    
    # aux2 = np.copy(aux)
    #aux.data[np.isnan(aux.data)] = -99999.0
    
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
    gdf['name'] = [ncObj[var].long_name]*len(data[3])
    units = ncObj[var].units
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in data[3]]
    
    if its != None:
        # t0 = ncObj['time'][0]
        t0 = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        # ti = ncObj['time'][its]
        ti = pd.Timedelta(seconds = int(ncObj['time'][its]))
        
        gdf['nTimeStep'] = [its]*len(data[3])
        # gdf['date'] = [str(ti.data).split('.')[0].replace('T', ' ')]*len(data[3])
        gdf['date'] = [str(t0 + ti)]*len(data[3])
        # gdf['nHours'] = [(ti.data - t0.data).item()/1e9]*len(data[3])
        gdf['nHours'] = [ti.total_seconds()/3600]*len(data[3])
        
    if subDom is not None:
        subDom = readSubDomain(subDom, epsg)
        gdf = gpd.clip(gdf, subDom)
    
    return gdf

def filledContours2gpd_mp(data, x, y, nv, name, units, longName, epoch, times, levels, epsg, subDom, its):
    ''' Run filledContours2gpd function with multiprocessing. NetCDF object can not be used with the Pool class,
        that's the reasson why only arrays are used as arguments.
        Parameters
            data: array
                values of the adcirc input file
            x: array
                longitude coordinates
            y: array
                latitude coordinates
            nv: array
                mesh elements
            name: string
                variable long name
            epoch: timestamp
                start date
            times: array
                time vector
            levels: np.array
                contour levels
            epsg: int
                coordinate system
            subDom: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            its: int
                timestep to be exported. If None (default) it is assumed that the file is not time-varying.
        Returns
            gdf: GeoDataFrame
                Linestrings as geometry and contour value in the "value" column.
                The name of the variable, the date and other info are also included in the columns
    '''
    nv = nv - 1 ## triangles starts from 1
    
    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    
    ## add max value over time and domain to the levels, just for polting porpuses.
    maxmax = np.max(data)
    if maxmax > levels[-1]:
        levels = np.append(levels, [np.ceil(maxmax)])
    
    aux = data[its, :]
    
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
    gdf['variable'] = [name]*len(data[3])
    gdf['name'] = [longName]*len(data[3])
    gdf['labelCol'] = [f'{x:0.2f} {units}' for x in data[3]]
    
    ti = pd.Timedelta(seconds = times[its])
        
    gdf['nTimeStep'] = [its]*len(data[3])
    gdf['date'] = [str(epoch + ti)]*len(data[3])
    gdf['nHours'] = [ti.total_seconds()/3600]*len(data[3])
        
    if subDom is not None:
        subDom = readSubDomain(subDom, epsg)
        gdf = gpd.clip(gdf, subDom)
    
    return gdf

def checkTimeVarying(ncObj):
    ''' Check if an adcirc input is time-varying or not.
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
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

def runExtractContours(ncObj, var, levels, conType, epsg, subDom, npro):
    ''' Run "contours2gpd" or "filledContours2gpd" if npro = 1 or "contours2gpd_mp" or "filledContours2gpd_mp" if npro > 1.
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
                Adcirc input file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            conType: string
                'polyline' or 'polygon'
            epsg: int
                coordinate system
            subDom: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
            npro: int
                number of worker processes. More info: https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool
        Returns
            gdf: GeoDataFrame
                Polygons or polylines as geometry columns. If the requested file is time-varying the GeoDataFrame will include all timesteps.
            
    '''
    timeVar = checkTimeVarying(ncObj)
    if timeVar == 0:
        if conType == 'polyline':
            gdf = contours2gpd(ncObj, var, levels, epsg, None, subDom)
        elif conType == 'polygon':
            gdf = filledContours2gpd(ncObj, var, levels, epsg, None, subDom)
        else:
            print('only "polyline" and "polygon" types are supported!')
            sys.exit(-1)
    else:
        if npro == 1:
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
        else: ## multiprocessing
            epoch = pd.to_datetime(ncObj['time'].units.split('since ')[1])
            fArgs = (ncObj[var][:, :], ncObj['x'][:], ncObj['y'][:], 
                    ncObj['element'][:,:], ncObj[var].name, ncObj[var].units,
                    ncObj[var].long_name, epoch, ncObj['time'][:], levels, epsg, subDom)
            
            fArgs = [fArgs + (it,) for it in range(ncObj['time'].size)]
            
            if npro == 999:
                p = Pool()
            else:
                p = Pool(npro)
            
            if conType == 'polyline':
                listGdf = p.starmap(contours2gpd_mp, tqdm(fArgs, total = len(fArgs)))
                p.close()
                p.join()
                
            elif conType == 'polygon':
                listGdf = p.starmap(filledContours2gpd_mp, tqdm(fArgs, total = len(fArgs)))
                p.close()
                p.join()
            else:
                print('only "polyline" and "polygon" types are supported!')
                sys.exit(-1)
            
        gdf = gpd.GeoDataFrame(pd.concat(listGdf, axis=0, ignore_index=True), crs=epsg)
        
    return gdf

def mesh2gdf(ncObj, epsg):
    ''' Write adcirc mesh as GeoDataFrame and extract centroid of each element. Used to create submesh
        Parameters:
            ncObj: netCDF4._netCDF4.Dataset
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
    ''' Generate a subdomain of the original adcirc input mesh based on a user input. NOT USED
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
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
    ''' Generate a subdomain of the original adcirc input mesh based on a subdomain. NOT USED
        Parameters
            ncObj: netCDF4._netCDF4.Dataset
                adcirc output
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
    ''' Write netcdf as xarray dataset. May be redundant but I had some problems reading adcirc files with xarray. NOT USED
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
    
def nc2shp(ncFile, var, levels, conType, epsg, pathOut, npro=1, subDomain=None, dateSubset=None):
    ''' Run all necesary functions to export adcirc outputs as shapefiles.
        Parameters
            ncFile: string
                path of the adcirc output, must be a netcdf file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            conType: string
                'polyline' or 'polygon'
            epsg: int
                coordinate system
            pathout: string
                complete path of the output file (*.shp or *.gpkg)
            npro: int
                number of worker processes. More info: https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool.
                For using all available processors, input 999.
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates
    '''
    # nc = nc2xr(ncFile, var)
    nc = netcdf.Dataset(ncFile, 'r')
    levels = np.arange(levels[0], levels[1], levels[2])

    gdf = runExtractContours(nc, var, levels, conType, epsg, subDomain, npro)
        
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
            fileName: string
                name of the colorbar img. Default tempColorbar.jpg.
            filePath: string
                path where the colorbar will be saved.
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
        screen1.icon.href = colorbarFile
        screen1.overlayxy = simplekml.OverlayXY(x= 0 , y = 0, xunits = simplekml.Units.fraction,
                                         yunits = simplekml.Units.fraction)
        screen1.screenxy = simplekml.ScreenXY(x = 0, y = 0.1, xunits = simplekml.Units.fraction,
                                         yunits = simplekml.Units.fraction)
        screen1.size.x = 0.08
        screen1.size.y = 0.55
        screen1.size.xunits = simplekml.Units.fraction
        screen1.size.yunits = simplekml.Units.fraction
    
    if logo == True:
        if logoFile == 'logo.png':
            aux0 = __file__
            aux1 = aux0.split('\\')
            aux2 = '\\'.join(aux1[:-2])
            logoFile = os.path.join(aux2, 'documentation', 'logo', logoFile)
            
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
    
    for i in gdf.index:
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

    for i in gdf.index:
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
    
def nc2kmz(ncFile, var, levels, conType, epsg, pathOut, npro=1, subDomain=None, overlay=True, logoFile='logo.png', colorbarFile='tempColorbar.jpg', 
            dateSubset=None, cmap='viridis'):
    ''' Run all necesary functions to export adcirc outputs as kmz.
        Parameters
            ncFile: string
                path of the adcirc output, must be a netcdf file
            var: string
                Name of the variable to export
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            conType: string
                'polyline' or 'polygon'
            epsg: int
                coordinate system
            pathout: string
                complete path of the output file (*.shp or *.gpkg)
            npro: int
                number of worker processes. More info: https://docs.python.org/3/library/multiprocessing.html#multiprocessing.pool.Pool
                For using all available processors, input 999.
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates                
            overlay: boolean. Default True
                If true overlay layer with logo and colorbar are added to the kmz file
            logoFile: string
                path of the logo image
            colorbarFile: string
                name of the colorbar img
            cmap: string
                name of the colormap
    '''
    # nc = nc2xr(ncFile, var)
    nc = netcdf.Dataset(ncFile, 'r')
    
    if checkTimeVarying(nc) == 1:
        print('Time-varying files can not be exported as kmz!')
        sys.exit(-1)
    else:
        levels = np.arange(levels[0], levels[1], levels[2])
        
        gdf = runExtractContours(nc, var, levels, conType, epsg, subDomain, npro)
        
        if conType == 'polygon':
            kml = polys2kml(gdf, levels, cmap)
        elif conType == 'polyline':
            kml = lines2kml(gdf, levels, cmap)
        else:
            print('Only "polygon" o "polyline" formats are supported')
            sys.exit(-1)
        
        if overlay == True:
            name = nc[var].long_name.capitalize()
            units = nc[var].units
            createColorbar(levels, name, units, cmap='viridis', fileName='tempColorbar.jpg', filePath='.')
            kmlScreenOverlays(kml, colorbar=True, colorbarFile='tempColorbar.jpg', logo=True, 
                        logoFile='logo.png', logoUnits='fraction', logoDims=None)
            
        kml.savekmz(pathOut, format = False)
        os.remove('tempColorbar.jpg')

