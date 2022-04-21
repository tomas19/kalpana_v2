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
        aux = ncObj[var][:].data#.to_masked_array()

    else:
        aux = ncObj[var][its, :].data#.to_masked_array()
    
    aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1)
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
    aux = data[its, :].data
    aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1)
        
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
        aux = ncObj[var][:].data#.to_masked_array()

    else:
        aux = ncObj[var][its, :].data#.to_masked_array()
    
    aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1)
    # aux2 = np.copy(aux)
    #aux.data[np.isnan(aux.data)] = -99999.0
    
    contoursf = plt.tricontourf(tri, aux, levels = levels)
    plt.close()  
    
    geoms = []
    
    ta = time.time()
    for icoll, coll in enumerate(contoursf.collections):
        ## get min an max malue of the interval
        vmin, vmax = contoursf.levels[icoll:icoll+2]
        vmean = np.mean([vmin, vmax])
        
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
    
    print(f'{(time.time() - ta)/60} min')
    
    ta = time.time()
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
        
    print(f'{(time.time() - ta)/60} min')
    
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
    
    aux = data[its, :].data
    aux = np.nan_to_num(aux, nan = -99999.0).reshape(-1)
    
    contoursf = plt.tricontourf(tri, aux, levels = levels)
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