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
import warnings
warnings.filterwarnings("ignore")
import time
import dask
from tqdm.dask import TqdmCallback
from vyperdatum.points import VyperPoints
from scipy.spatial import KDTree

def gdfChangeVerUnit(gdf, ini, out):
    ''' Change vertical units of the columns of a GeoDataFrame which
        contains "z" on the name
        Parameters
            gdf: GeoDataFrame
                out of the runExtractContours function
            ini: string
                initial unit eg: 'm' or 'ft'
            out: string
                final unit eg: 'm' or 'ft'
        Returns
            gdf2: GeoDataFrame
                updated geodataframe with the new unit added on the name 
                to the columns with "z" in it.
    '''
    gdf2 = gdf.copy()
    colsIni = [x if 'z' not in x else f'{x}_{ini}' for x in gdf2.columns]
    gdf2.columns = colsIni
    
    cols = [x for x in gdf2.columns if 'z' in x]
    if ini == 'm' and out == 'ft':
        for col in cols:
            # colnew = col.replace(ini, out)
            # gdf2[colnew] = gdf2[col] * 3.2808399
            gdf2[col] = gdf2[col] * 3.2808399
    elif ini == 'ft' and out == 'm':
        for col in cols:
            # colnew = col.replace(ini, out)
            # gdf2[colnew] = gdf2[col] / 3.2808399
            gdf2[col] = gdf2[col] / 3.2808399
            
    return gdf2

def dzDatums(vdatum_directory, x, y, epsg, vdatumIn, vdatumOut):
    ''' Get the vertical difference of two datums on a group of locations
        Parameters
            vdatum_directory: string
                full path of the instalation folder of vdatum (https://vdatum.noaa.gov/)
            x, y: numpy array
                x and y-coordinate of the group of points
            epsg: int
                coordinate system code of the data. eg: 4326 for wgs84 (lat/lon)
            vdatumIn, vdatumOut: string
                name of the input and output vertical datums. Mean sea level is "tss"
                For checking the available datums:
                from vyperdatum.pipeline import datum_definition
                list(datum_definition.keys())
        Returns
            df: dataframe
                dataframe with the coordinates of the requested points. The dz0 column is the raw
                datum difference from vdatum, but it can also output NaN values. In this case the difference
                of the closest valid point is assigned used the scipy KDTree algorithm.
                
    '''
    vp = VyperPoints(vdatum_directory = vdatum_directory)
    z = np.zeros(len(x))
    vp = VyperPoints(silent = True)
    vp.transform_points((epsg, vdatumIn), (epsg, vdatumOut), x, y, z = z)
    df = pd.DataFrame({'x': x, 'y': y, 'dz0': vp.z})
    indexVal = np.where(~np.isnan(df['dz0'].values))[0]
    indexNan = np.where(np.isnan(df['dz0'].values))[0]
    tree = KDTree(list(zip(df.loc[indexVal, 'x'], df.loc[indexVal, 'y'])))
    query = tree.query(list(zip(df.loc[indexNan, 'x'], df.loc[indexNan, 'y'])))
    aux0 = df.loc[indexVal, 'dz0']
    aux1 = aux0.values[query[1]]
    df['dz1'] = [np.nan]*len(df)
    df.loc[indexNan, 'dz1'] = aux1
    df['dz'] = df['dz0'].fillna(0) + df['dz1'].fillna(0)
    
    return df
    
def gdfChangeVertDatum(vdatum_directory, gdf, vdatumIn = 'tss', vdatumOut = 'navd88'):
        ''' Change the vertical datum of the columns of a GeoDataFrame which
            contains "z" on the name
            Parameters
                vdatum_directory: string
                    full path of the instalation folder of vdatum (https://vdatum.noaa.gov/)
            gdf: GeoDataFrame
                out of the runExtractContours function
            vdatumIn, vdatumOut: string. Default 'tss' and 'navd88'
                name of the input and output vertical datums. Mean sea level is "tss"
                For checking the available datums:
                from vyperdatum.pipeline import datum_definition
                list(datum_definition.keys())
            gdf2: GeoDataFrame
                updated geodataframe with the new datum added on the name 
                to the columns with "z" on it.
        '''
        ## get points
        gdf2 = gdf.copy()
        dummyList = [list(gdf2['geometry'][x].centroid.coords)[0] for x in gdf2.index]
        aux = list(zip(*dummyList))
        x = np.array(aux[0])
        y = np.array(aux[1])
        z = np.zeros(len(x))
        epsg = int(gdf2.crs.to_string().split(':')[1])
        
        # gdf2['xPolyCentroid'] = x
        # gdf2['yPolyCentroid'] = y
        
        ## use vyperdatum to get difference betweem vdatumIn and vdatumOut
        ## dz in m
        df = dzDatums(vdatum_directory, x, y, epsg, vdatumIn, vdatumOut)
        
        newCols = [x if 'z' not in x else f'{x}_{vdatumIn}' for x in gdf2.columns]
        gdf2.columns = newCols
        aux = [x for x in gdf2.columns if vdatumIn in x]
        for col in aux:
            # newcol = col.replace(vdatumIn, vdatumOut)
            # gdf2[newcol] = gdf2[col] + df.dz
            gdf2[col] = gdf2[col] + df.dz
        
        return gdf2

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
    areas = np.array([signedArea(p) for p in polys])
    outer = [polys[x] for x in list(np.where(areas >= 0)[0])]
    inner = [polys[x] for x in list(np.where(areas < 0)[0])]
    
    return outer, inner

def signedArea(ring):
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
    sArea = np.cross(ring, v2).sum() / 2.0
    
    return sArea

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

def filledContours2gpd(tri, data, levels, epsg):
    ''' Dataset to GeoDataFrame with filled contours as shapely polygons.
        A dask delayed python decorator is used to handle the code's parallelization.
        Parameters
            tri: mpl.tri.Triangulation
                triangulation of the adcirc mesh
            data: numpy array
                1D array with the data of the adcirc output variable
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            epsg: int
                coordinate system
        Returns
            gdf: GeoDataFrame
                Polygons as geometry and contours min, average, max values as columns. 
                The name of the variable, the date and other info are also included in the columns
    '''
        
    contoursf = plt.tricontourf(tri, data, levels = levels)
    plt.close()  
    
    @dask.delayed
    def getGeom(icoll, coll):
        vmin, vmax = contoursf.levels[icoll:icoll+2]
        vmean = np.mean([vmin, vmax])  
        geoms = []
        
        for p in coll.get_paths():
            p.simplify_threshold = 0.0
            # Removing polygons with less than 3 vertices
            polys = [g for g in p.to_polygons() if g.shape[0] >= 3] 
            # classify polygons based on the signed area
            outer, inner = classifyPolygons(polys)

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
        
        return geoms
    
    tasks = [getGeom(icoll, coll) for icoll, coll in enumerate(contoursf.collections)]
    with TqdmCallback(desc = "Compute contours using Dask"):
        daskGeoms = dask.compute(tasks, scheduler = 'threads')
    
    geoms = list(itertools.chain(*daskGeoms[0]))
    
    ## define geopandas
    data = list(zip(*geoms))
    gdf = gpd.GeoDataFrame(crs = epsg, geometry = list(data[0]), 
                           index = range(len(data[0])))
    gdf['zMin'] = data[1]
    gdf['zMean'] = data[3]
    gdf['zMax'] = data[2]
    
    return gdf

def contours2gpd(tri, data, levels, epsg):
    ''' Dataset to GeoDataFrame with contours as shapely LineStrings
        A dask delayed python decorator is used to handle the code's parallelization.
        Parameters
            tri: mpl.tri.Triangulation
                triangulation of the adcirc mesh
            data: numpy array
                1D array with the data of the adcirc output variable
            levels: np.array
                Contour levels. The max value in the entire doman and over all timesteps is added to the requested levels.
            epsg: int
                coordinate system
        Returns
            gdf: GeoDataFrame
                Linestrings as geometry and contour value in the "value" column.
                The name of the variable, the date and other info are also included in the columns
    '''

    contours = plt.tricontour(tri, data, levels = levels)
    plt.close()
    
    @dask.delayed
    def getGeom(icoll, coll):

        ## iteration over lines
        geoms = []

        for icon, con in enumerate(contours.collections):
            val = contours.levels[icon]
            paths = con.get_paths()
            ## dismiss lines with less than 2 vertices
            aux0 = [(LineString(path.vertices), val) for path in paths if len(path.vertices) > 2]
            geoms.append(aux0)
        return geoms

    tasks = [getGeom(icoll, coll) for icoll, coll in enumerate(contours.collections)]
    with TqdmCallback(desc = "Compute contours using Dask"):
        daskGeoms = dask.compute(tasks, scheduler = 'threads')

    geoms = list(itertools.chain(*daskGeoms[0]))
    geoms = list(itertools.chain(*geoms))
    
    ## define geopandas
    data = list(zip(*geoms))
    gdf = gpd.GeoDataFrame(crs = epsg, geometry = list(data[0]), 
                           index = range(len(data[0])))
   
    gdf['z'] = data[1]
    
    return gdf

def runExtractContours(ncObj, var, levels, conType, epsg):
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
        Returns
            gdf: GeoDataFrame
                Polygons or polylines as geometry columns. If the requested file is time-varying the GeoDataFrame will include all timesteps.
            
    '''
    
    nv = ncObj['element'][:,:] - 1 ## triangles starts from 1
    x = ncObj['x'][:].data
    y = ncObj['y'][:].data
    
    vname = ncObj[var].name
    lname = ncObj[var].long_name
    u = ncObj[var].units

    ## matplotlib triangulation
    tri = mpl.tri.Triangulation(x, y, nv)
    
    timeVar = checkTimeVarying(ncObj)
    
    if timeVar == 0:
        aux = ncObj[var][:].data
        aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1) 
        
        if conType == 'polyline':
            labelCol = 'z'
            gdf = contours2gpd(tri, aux, levels, epsg)
        
        elif conType == 'polygon':
            labelCol = 'zMean'
            ## add max value over time and domain to the levels, just for polting porpuses.
            maxmax = np.max(ncObj[var])
            if maxmax > levels[-1]:
                levels = np.append(levels, [np.ceil(maxmax)])
            
            gdf = filledContours2gpd(tri, aux, levels, epsg)
        
        else:
            print('only "polyline" and "polygon" types are supported!')
            sys.exit(-1)
            
        gdf['variable'] = [vname]*len(gdf)
        gdf['name'] = [lname]*len(gdf)
        gdf['labelCol'] = [f'{x:0.2f} {u}' for x in gdf[labelCol]]
            
    else:
        t0 = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        listGdf = []
        if conType == 'polyline':
            labelCol = 'z'
            for t in tqdm(range(ncObj['time'].shape[0])):
                aux = ncObj[var][t, :].data
                aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1) 
                gdfi = contours2gpd(tri, aux, levels, epsg)
                
                ti = pd.Timedelta(seconds = int(ncObj['time'][t]))
                gdfi['nTimeStep'] = [t]*len(gdfi)
                gdfi['date'] = [str(t0 + ti)]*len(gdfi)
                gdfi['nHours'] = [ti.total_seconds()/3600]*len(gdfi)
                listGdf.append(gdfi)
        
        elif conType == 'polygon':
            labelCol = 'zMean'
            for t in tqdm(range(ncObj['time'].shape[0])):
                aux = ncObj[var][t, :].data
                aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1) 
                gdfi = filledContours2gpd(tri, aux, levels, epsg)
                
                ti = pd.Timedelta(seconds = int(ncObj['time'][t]))
                gdfi['nTimeStep'] = [t]*len(gdfi)
                gdfi['date'] = [str(t0 + ti)]*len(gdfi)
                gdfi['nHours'] = [ti.total_seconds()/3600]*len(gdfi)
                listGdf.append(gdfi)
        else:
            print('only "polyline" and "polygon" types are supported!')
            sys.exit(-1)
            
        gdf = gpd.GeoDataFrame(pd.concat(listGdf, axis=0, ignore_index=True), crs=epsg)
        gdf['variable'] = [vname]*len(gdf)
        gdf['name'] = [lname]*len(gdf)
        gdf['labelCol'] = [f'{x:0.2f} {u}' for x in gdf[labelCol]]
        
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
    gdf = gpd.GeoDataFrame(geometry = pols, crs = epsg)
    gdf['centX'] = xvertices.mean(axis = 1)
    gdf['centY'] = yvertices.mean(axis = 1)
    gdf['v1'] = nv.data[:, 0]
    gdf['v2'] = nv.data[:, 1]
    gdf['v3'] = nv.data[:, 2]
    
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
        if subDomain.endswith('.shp') or subDomain.endswith('.gpkg'):
            gdfSubDomain = gpd.read_file(subDomain)
        elif subDomain.endswith('.kml'):
            gpd.io.file.fiona.drvsupport.supported_drivers['KML'] = 'rw'
            gdfSubDomain = gpd.read_file(subDomain, driver = 'KML')
        else:
            print('Only shape, geopackage and kml formats are suported for sub domain generation!')
            sys.exit(-1)
        xAux, yAux = gdfSubDomain.geometry[0].exterior.coords.xy
        extCoords = list(zip(xAux, yAux))
    
    elif type(subDomain) == list and len(subDomain) == 4:
        ## only UL lon, UL lat, LR lon LR lat
        ulLon, ulLat, lrLon, lrLat = subDomain
        extCoords = [(ulLon, ulLat), (ulLon, lrLat), (lrLon, lrLat), (lrLon, ulLat), (ulLon, ulLat)]
        poly = Polygon(extCoords)
        gdfSubDomain = gpd.GeoDataFrame(geometry = [poly], crs = epsg)
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
        
        ds[var].attrs['long_name'] = ncObj[var].long_name.capitalize()
        ds[var].attrs['units'] = ncObj[var].units
        # ds['time'].attrs['base_date'] = pd.to_datetime(ncObj['time'].units.split('since ')[1])
        # ds['time'].attrs['units'] = f'Seconds since {epoch}'
        
    return ds
    
def nc2shp(ncFile, var, levels, conType, epsgIn, epsgOut, vUnitIn, vUnitOut, vDatumIn, vDatumOut, pathOut, vDatumPath,
           subDomain=None):
    ''' Run all necesary functions to export adcirc outputs as shapefiles.
        Parameters
            ncFile: string
                path of the adcirc output, must be a netcdf file
            var: string
                Name of the variable to export
            levels:list
                Contour levels. Min, Max and Step. Max is not included as in np.arange method.
            conType: string
                'polyline' or 'polygon'
            epsgIn: int
                coordinate system of the adcirc input
            epsgOut: int
                coordinate system of the output shapefile
            vUnitIn, vUnitOut: string
                input and output vertical units. For the momment only supported 'm' and 'ft'
            vdatumIn, vdatumOut: string
                name of the input and output vertical datums. Mean sea level is "tss"
                For checking the available datums:
                from vyperdatum.pipeline import datum_definition
                list(datum_definition.keys())
            pathout: string
                complete path of the output file (*.shp or *.gpkg)
            vdatum_directory: string
                full path of the instalation folder of vdatum (https://vdatum.noaa.gov/)
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates. The crs must be the same of the
                adcirc input file.
    '''
    
    nc = netcdf.Dataset(ncFile, 'r')
    levels = np.arange(levels[0], levels[1], levels[2])

    gdf = runExtractContours(nc, var, levels, conType, epsgIn)
    
    if subDomain is not None:
        subDom = readSubDomain(subDomain, epsgIn)
        gdf = gpd.clip(gdf, subDom)

    if vDatumIn == vDatumOut:
        pass
    else:
        gdf = gdfChangeVertDatum(vDatumPath, gdf, vDatumIn, vDatumOut)
    
    if vUnitIn == vUnitOut:
        pass
    else:
        gdf = gdfChangeVerUnit(gdf, vUnitIn, vUnitOut)
    
    if epsgIn == epsgOut:
        pass
    else:
        gdf = gdf.to_crs(epsgOut)
        
    if pathOut.endswith('.shp'):
        gdf.to_file(pathOut)
    elif pathOut.endswith('.gpkg'):
        gdf.to_file(pathOut, driver = 'GPKG')
       
    elif pathOut.endswith('.wkt'):
        gdf.to_csv(pathOut)
    
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
            logoFile = os.path.join(aux2, 'documentation', 'logoForKmz', logoFile)
            
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
        value = gdf.loc[i, 'z']
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
            value = gdf.loc[i, 'zMean']
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
                value = gdf.loc[i, 'zMean']
                r, g, b, a = m.to_rgba(value)
                col = simplekml.Color.rgb(int(255*r), int(255*g), int(255*b))
                pol.style.linestyle.color = col
                pol.style.linestyle.width = 2
                pol.description = gdf.loc[i, 'labelCol']
        #         pol.style.polystyle.color = col
                pol.style.polystyle.color = simplekml.Color.changealphaint(100, col)

    return kml

def countVertices(gdf):
    ''' Count vertices of the polygons of a GeoDataFrame.
        Parameters:
            gdf: GeoDataFrame
                output of runExtractContours function
        Returns:
            nVertices: list
                number of vertices of each polygon
    '''
    nVertices = []
    for ix, x in enumerate(gdf.geometry):
        poly = x
        a = 0
        try:
            a = a + len(list(x.exterior.coords))
        except:
            for pol in poly:
                a = a + len(list(pol.exterior.coords))
        nVertices.append(a)
        
        return nVertices

def katana(geometry, threshold, count=0):
    ''' Split a Polygon into two or more parts across it's shortest dimension
        function taken from: 
        https://snorfalorpagus.net/blog/2016/03/13/splitting-large-polygons-for-faster-intersections
        Parameters:
            geometry: shapely polygon
                polygon to devide
            threshold: float
                maximum area of the subpolygons
            count: int
                number of subareas. Default 0.
        Returns
            final_result: list
                list with new polygons
    '''
    bounds = geometry.bounds
    width = bounds[2] - bounds[0]
    height = bounds[3] - bounds[1]
    if max(width, height) <= threshold or count == 250:
        # either the polygon is smaller than the threshold, or the maximum
        # number of recursions has been reached
        return [geometry]
    if height >= width:
        # split left to right
        a = box(bounds[0], bounds[1], bounds[2], bounds[1]+height/2)
        b = box(bounds[0], bounds[1]+height/2, bounds[2], bounds[3])
    else:
        # split top to bottom
        a = box(bounds[0], bounds[1], bounds[0]+width/2, bounds[3])
        b = box(bounds[0]+width/2, bounds[1], bounds[2], bounds[3])
    result = []
    for d in (a, b,):
        c = geometry.intersection(d)
        if not isinstance(c, GeometryCollection):
            c = [c]
        for e in c:
            if isinstance(e, (Polygon, MultiPolygon)):
                result.extend(katana(e, threshold, count+1))
    if count > 0:
        return result
    # convert multipart into singlepart
    final_result = []
    for g in result:
        if isinstance(g, MultiPolygon):
            final_result.extend(g)
        else:
            final_result.append(g)
    
    return final_result
    
def nc2kmz(ncFile, var, levels, conType, epsg, pathOut,subDomain=None, overlay=True, logoFile='logo.png', colorbarFile='tempColorbar.jpg', 
           cmap='viridis', thresVertices=20_000):
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
            subDomain: str or list. Default None
                complete path of the subdomain polygon kml or shapelfile, or list with the
                uper-left x, upper-left y, lower-right x and lower-right y coordinates The crs must be the same of the
                adcirc input file.
            overlay: boolean. Default True
                If true overlay layer with logo and colorbar are added to the kmz file
            logoFile: string. Default 'logo.png'
                path of the logo image
            colorbarFile: string. Default 'tempColorbar.jpg'
                name of the colorbar img
            cmap: string. Default 'viridis'
                name of the colormap
            thresVertices: int
                maximum number of vertices allowed per polygon. If a polygon has more vertices, the
                katana function will be used.
    '''
    # nc = nc2xr(ncFile, var)
    nc = netcdf.Dataset(ncFile, 'r')
    
    if checkTimeVarying(nc) == 1:
        print('Time-varying files can not be exported as kmz!')
        sys.exit(-1)
    else:
        levels = np.arange(levels[0], levels[1], levels[2])
        
        gdf = runExtractContours(nc, var, levels, conType, epsg)
        
        # if subDomain is not None:
            # subDom = readSubDomain(subDomain, epsg)
            # gdf = gpd.clip(gdf, subDom)
        
        # if conType == 'polygon':
            # kml = polys2kml(gdf, levels, cmap)
        # elif conType == 'polyline':
            # kml = lines2kml(gdf, levels, cmap)
        # else:
            # print('Only "polygon" o "polyline" formats are supported')
            # sys.exit(-1)
        
        # if overlay == True:
            # name = nc[var].long_name.capitalize()
            # units = nc[var].units
            # createColorbar(levels, name, units, cmap='viridis', fileName='tempColorbar.jpg', filePath='.')
            # kmlScreenOverlays(kml, colorbar=True, colorbarFile='tempColorbar.jpg', logo=True, 
                        # logoFile='logo.png', logoUnits='fraction', logoDims=None)
            
        # kml.savekmz(pathOut, format = False)
        # os.remove('tempColorbar.jpg')
        
        return gdf

