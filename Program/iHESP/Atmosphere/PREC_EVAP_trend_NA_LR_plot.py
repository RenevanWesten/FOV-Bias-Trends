#Program plot the precipitation minus evaporation trend

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
directory	            = '../../../Data/iHESP/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

experiment	= 'PI_control_LR'
#experiment	= 'RCP_LR'

if experiment == 'PI_control_LR':
	fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_NA_trend_year_21-100_PI_control_LR.nc', 'r')
if experiment == 'RCP_LR':
	fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_NA_trend_year_2000-2100_RCP_LR.nc', 'r')

lon_atm		= fh.variables['lon'][:] 	
lat_atm		= fh.variables['lat'][:] 
PREC_trend	= fh.variables['PREC_trend'][:] 	
EVAP_trend	= fh.variables['EVAP_trend'][:] 					
P_E_trend	= fh.variables['P_E_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_atm, lat_atm, P_E_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('P-E trend (mm day$^{-1}$ per century)')

gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-60, -30])
gl1.ylocator = mticker.FixedLocator([20, 30, 40])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False)
ax.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_LR':
	ax.set_title('b) P-E trend, LR-CESM, PI control (21 - 100)')

if experiment == 'RCP_LR':
	ax.set_title('P-E trend, LR-CESM, Hist/RCP8.5 (2000 - 2100)')
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_atm, lat_atm, PREC_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Precipitation trend (mm day$^{-1}$ per century)')

gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-60, -30])
gl1.ylocator = mticker.FixedLocator([20, 30, 40])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False)
ax.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_LR':
	ax.set_title('d) Precipitation trend, LR-CESM, PI control (21 - 100)')

if experiment == 'RCP_LR':
	ax.set_title('Precipitation trend, LR-CESM, Hist/RCP8.5 (2000 - 2100)')
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_atm, lat_atm, EVAP_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Evaporation trend (mm day$^{-1}$ per century)')

gl1	     = ax.gridlines(draw_labels=True, dms = True, x_inline=False, y_inline=False, linewidth = 0.0)
gl1.top_labels = False
gl1.right_labels = False
gl1.xlocator = mticker.FixedLocator([-60, -30])
gl1.ylocator = mticker.FixedLocator([20, 30, 40])
gl1.xlabel_style = {'rotation':0}
gl2 	= ax.gridlines(draw_labels=False, dms = True, x_inline=False, y_inline=False)
ax.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_LR':
	ax.set_title('f) Evaporation trend, LR-CESM, PI control (21 - 100)')

if experiment == 'RCP_LR':
	ax.set_title('Evaporation trend, LR-CESM, Hist/RCP8.5 (2000 - 2100)')

show()


