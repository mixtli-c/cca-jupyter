import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature

fname = 'F:\\jupyter\\shapefiles\\st_mx.shp'

ax = plt.axes(projection=ccrs.Robinson())
shape_feature = ShapelyFeature(Reader(fname).geometries(),
                                ccrs.PlateCarree(), facecolor='none', edgecolor='black')
ax.set_extent((-105,-95,16,22))

ax.add_feature(shape_feature)
plt.show()
