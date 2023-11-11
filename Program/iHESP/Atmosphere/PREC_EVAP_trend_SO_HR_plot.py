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

#Making pathway to folder with all data
directory	            = '../../../Data/iHESP/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

experiment	= 'PI_control_HR'
#experiment	= 'RCP_HR'


if experiment == 'PI_control_HR':
	fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_SO_trend_year_21-100_PI_control_HR.nc', 'r')
if experiment == 'RCP_HR':
	fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_SO_trend_year_2000-2100_RCP_HR.nc', 'r')


lon_atm		= fh.variables['lon'][:] 	
lat_atm		= fh.variables['lat'][:] 
PREC_trend	= fh.variables['PREC_trend'][:] 	
EVAP_trend	= fh.variables['EVAP_trend'][:] 					
P_E_trend	= fh.variables['P_E_trend'][:]
P_E_trend_sig	= fh.variables['P_E_trend_sig'][:] 	

fh.close()

P_E_trend_sig	= ma.masked_where(P_E_trend_sig < 0.95, P_E_trend_sig)
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_atm, lat_atm, P_E_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

CS2 	= ax.contourf(lon_atm, lat_atm, P_E_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '\\', transform=ccrs.PlateCarree())
CS2 	= ax.contourf(lon_atm, lat_atm, P_E_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '//', transform=ccrs.PlateCarree())	

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('P-E trend (mm day$^{-1}$ per century)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_HR':
	ax.set_title('a) P-E trend, HR-CESM, PI control (21 - 100)')

if experiment == 'RCP_HR':
	ax.set_title('c) P-E trend, HR-CESM, Hist/RCP8.5 (2000 - 2100)')

#-----------------------------------------------------------------------------------------


fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_atm, lat_atm, PREC_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Precipitation trend (mm day$^{-1}$ per century)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_HR':
	ax.set_title('c) Precipitation trend, HR-CESM, PI control (21 - 100)')

if experiment == 'RCP_HR':
	ax.set_title('Precipitation trend, HR-CESM, Hist/RCP8.5 (2000 - 2100)')
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_atm, lat_atm, EVAP_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Evaporation trend (mm day$^{-1}$ per century)')

gl = ax.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
ax.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax.coastlines('50m')
ax.add_feature(cfeature.LAND, zorder=0)

if experiment == 'PI_control_HR':
	ax.set_title('e) Evaporation trend, HR-CESM, PI control (21 - 100)')

if experiment == 'RCP_HR':
	ax.set_title('Evaporation trend, HR-CESM, Hist/RCP8.5 (2000 - 2100)')

show()


