#Program the PI control LR-CESM trends over the entire first 100 years

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
directory	= '../../../Data/iHESP/'

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_SO_trend_year_1-100_PI_control_LR.nc', 'r')

lon_atm		= fh.variables['lon'][:] 	
lat_atm		= fh.variables['lat'][:] 					
P_E_trend	= fh.variables['P_E_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_atm, lat_atm, P_E_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

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
ax.set_title('a) P-E trend, LR-CESM, PI control (1 - 100)')

#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100

fh = netcdf.Dataset(directory+'Ocean/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_1-100_PI_control_LR.nc', 'r')

lon_ocn		= fh.variables['TLONG'][:] 	
lat_ocn		= fh.variables['TLAT'][:] 		
salt_trend	= fh.variables['SALT_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})

CS      = ax.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

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
ax.set_title('b) Salinity (0 - 100 m) trend, LR-CESM, PI control (1 - 100)')

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/SST_NA_trend_year_1-100_PI_control_LR.nc', 'r')

lon_ocn		= fh.variables['TLONG'][:] 	
lat_ocn		= fh.variables['TLAT'][:] 		
temp_trend	= fh.variables['TEMP_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_ocn, lat_ocn, temp_trend, levels = np.arange(-2, 2.1, 0.1), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-2, 2.1, 1), cax=ax_cb)
cbar.set_label('SST trend ($^{\circ}$C per century)')

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
ax.set_title('c) SST trend, LR-CESM, PI control (1 - 100)')

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_NA_trend_year_1-100_PI_control_LR.nc', 'r')

lon_atm		= fh.variables['lon'][:] 	
lat_atm		= fh.variables['lat'][:] 					
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
ax.set_title('d) P-E trend, LR-CESM, PI control (1 - 100)')

#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100

fh = netcdf.Dataset(directory+'Ocean/SALT_NA_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_1-100_PI_control_LR.nc', 'r')

lon_ocn		= fh.variables['TLONG'][:] 	
lat_ocn		= fh.variables['TLAT'][:] 		
salt_trend	= fh.variables['SALT_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)')

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
ax.set_title('e) Salinity (0 - 100 m) trend, LR-CESM, PI control (1 - 100)')

#-----------------------------------------------------------------------------------------

depth_min	= 1000
depth_max	= 3000

fh = netcdf.Dataset(directory+'Ocean/SALT_NA_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_1-100_PI_control_LR.nc', 'r')

lon_ocn		= fh.variables['TLONG'][:] 	
lat_ocn		= fh.variables['TLAT'][:] 		
salt_trend	= fh.variables['SALT_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.Stereographic(central_longitude=-40, central_latitude=40)})

CS      = ax.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

divider = make_axes_locatable(ax)
ax_cb   = divider.new_horizontal(size="5%", pad=0.1, axes_class=plt.Axes)
fig.add_axes(ax_cb)

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), cax=ax_cb)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)')

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
ax.set_title('f) Salinity (1000 - 3000 m) trend, LR-CESM, PI control (1 - 100)')
#-----------------------------------------------------------------------------------------

show()
