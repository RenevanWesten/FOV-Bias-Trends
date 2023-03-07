#Program plots the salinity trends

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker

#Making pathway to folder with all data
directory	   = '../../../Data/iHESP/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min	= 0
depth_max	= 100

#-----------------------------------------------------------------------------------------	

fh = netcdf.Dataset(directory+'Ocean/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_2000-2100_RCP_HR.nc', 'r')

lon_ocn		    = fh.variables['lon'][:] 	
lat_ocn		    = fh.variables['lat'][:] 		
salt_trend	    = fh.variables['SALT_trend'][:] 
salt_trend_sig	= fh.variables['SALT_trend_sig'][:]	

fh.close()

salt_trend_sig	= ma.masked_where(salt_trend_sig < 0.95, salt_trend_sig)
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

CS2 	= ax.contourf(lon_ocn, lat_ocn, salt_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '\\', transform=ccrs.PlateCarree())
CS2 	= ax.contourf(lon_ocn, lat_ocn, salt_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '//', transform=ccrs.PlateCarree())	

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.01, 1), cax=ax_cb)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

ax.set_title('a) Salinity (0 - 100 m) trend, HR-CESM, Hist/RCP8.5 (2000 - 2100)')

show()

