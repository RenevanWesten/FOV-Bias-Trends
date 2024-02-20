#Program plots the surface (0-100m) salinity in the Indian Ocean

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats
from scipy.stats import genextreme
from matplotlib.colors import LogNorm
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.tri as tri

#Making pathway to folder with all data
directory       = '../../../Data/POP/'
directory_rean  = '../../../Data/Reanalysis/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100

#-----------------------------------------------------------------------------------------

#Read in the stand-alone POP
fh = netcdf.Dataset(directory+'Ocean/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m.nc', 'r')

lon		= fh.variables['lon'][:]	
lat		= fh.variables['lat'][:]	
salt_all= fh.variables['SALT'][:]		#Salinity

fh.close()

salt_all	= np.mean(salt_all, axis = 0)

#-----------------------------------------------------------------------------------------

#Read in the salinity field for the Reanalysis product
files	= glob.glob(directory_rean+'Data/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m/Reanalysis_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

for file_i in range(len(files)):
	print(files[file_i])
	fh = netcdf.Dataset(files[file_i], 'r')

	lon_rean	= fh.variables['lon'][:]	
	lat_rean	= fh.variables['lat'][:]	
	salt        = fh.variables['SALT'][:]		#Salinity

	fh.close()

	if file_i == 0:
		salt_rean	= ma.masked_all((len(files), len(lat_rean), len(lon_rean)))

	salt_rean[file_i]	= salt

salt_rean	= np.mean(salt_rean, axis = 0)
#-----------------------------------------------------------------------------------------#-----------------------------------------------------------------------------------------

mask_rean	       = salt_rean.mask == False
lon_rean, lat_rean = np.meshgrid(lon_rean, lat_rean)
lon_rean		   = lon_rean[mask_rean]
lat_rean		   = lat_rean[mask_rean]
salt_rean	       = salt_rean[mask_rean]

#Triangulate the grid of the reanalysis
triang            = tri.Triangulation(lon_rean, lat_rean)
salt_rean         = tri.LinearTriInterpolator(triang, salt_rean)
lon, lat          = np.meshgrid(lon, lat)
salt_rean         = salt_rean(lon, lat)
salt_rean         = ma.masked_array(salt_rean, mask = salt_all.mask)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize = (6.8, 4.8))

CS      = ax.contourf(lon, lat, salt_all, levels = np.arange(33, 37.1, 0.1), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(33, 37.1, 1), cax=ax_cb)
cbar.set_label('Salinity (g kg$^{-1}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

ax.set_title('d) Salinity (0 - 100 m), Stand-alone POP (245 - 274)')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Salinity anomaly plot
fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()}, figsize = (6.8, 4.8))

CS      = ax.contourf(lon, lat, salt_all - salt_rean, levels = np.arange(-1, 1.1, 0.1), extend = 'both', cmap = 'PuOr_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 0.5), cax=ax_cb)
cbar.set_label('Salinity difference (g kg$^{-1}$)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

ax.set_title('d) Salinity (0 - 100 m), Stand-alone POP (245 - 274)')

show()

