kalpana modules
1) Functions to export adcirc as shapefiles or kmz files.
2) Downscale flooding outputs aiming to increase the resolution near coastline. WIP

Usage: follow examples on the examples folder

TODO:
1) Fix warning GeoDataFrane index
2) Implement date subset
3) Avoid using geopandas clip to select subset. Maybe it is faster to select elements inside subdomain and then do the calculations
4) Wind data, magnitud and direction needs to be calculated from vectors
5) Impelement discrete colorbar for kmz
6) Add different colormap for topo
7) Look for a faster way to save shapefiles
