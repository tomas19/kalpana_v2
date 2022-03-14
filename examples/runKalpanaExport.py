import sys
import time
import warnings

warnings.filterwarnings("ignore")

now = time.time()
inpf = sys.argv[1]

with open(inpf, 'r') as inp:
    inputs = inp.readlines()[-1].split()

scriptPath = inputs[0]
ncFile = inputs[1]
var = inputs[2]
levels = [float(x) for x in inputs[3].split(',')]
conType = inputs[4]
epsg = int(inputs[5])
pathOut = inputs[6]

try:
    subDom = inputs[7]
except:
    subDom = 'NoSubDomain'

sys.path.append(scriptPath)
try:
    from kalpanaExport import nc2shp, nc2kmz
except:
    print('kalpana module not found!')
    sys.exit(-1)

if pathOut.endswith('.kmz'):
    if subDom != 'NoSubDomain':
        gdf = nc2kmz(ncFile, var, levels, conType, epsg, pathOut, subDomain=subDom)
    else:
        gdf = nc2kmz(ncFile, var, levels, conType, epsg, pathOut)

elif pathOut.endswith('.shp') or pathOut.endswith('.gpkp'):
    if subDom != 'NoSubDomain':
        gdf = nc2shp(ncFile, var, levels, conType, epsg, pathOut, subDomain=subDom)
    else:
        gdf = nc2shp(ncFile, var, levels, conType, epsg, pathOut)

else:
    print('Only ".kmz", ".shp" or ".gpkg" formats are suported!')
    sys.exit(-1)

print(f'Script finished succsesfully after: {(time.time() - now)/60:0.3f} min')
