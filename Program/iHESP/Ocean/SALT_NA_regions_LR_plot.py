#Program plots the depth-averaged salinity over the three basins

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
directory	   = '../../../Data/iHESP/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	salt_lab	= fh.variables['SALT_Lab'][:]	
	salt_irm	= fh.variables['SALT_Irm'][:]	
	salt_ice	= fh.variables['SALT_Ice'][:]	

	fh.close()

	return time, salt_lab, salt_irm, salt_ice

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, salt_lab_control, salt_irm_control, salt_ice_control	= ReadinData(directory+'Ocean/SALT_NA_regions_PI_control_LR.nc')
time_rcp, salt_lab_rcp, salt_irm_rcp, salt_ice_rcp			= ReadinData(directory+'Ocean/SALT_NA_regions_RCP_LR.nc')
time_control								+= 1599
#-----------------------------------------------------------------------------------------

depth_min	= 0
depth_max	= 100

fh = netcdf.Dataset(directory+'Ocean/SALT_NA_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_21-100_PI_control_LR.nc', 'r')

lon_ocn			= fh.variables['TLONG'][:] 	
lat_ocn			= fh.variables['TLAT'][:] 		
salt_0_100_trend= fh.variables['SALT_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

depth_min	= 1000
depth_max	= 3000

fh = netcdf.Dataset(directory+'Ocean/SALT_NA_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_21-100_PI_control_LR.nc', 'r')
		
salt_1000_3000_trend	= fh.variables['SALT_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_lab	= plot(time_control, salt_lab_control, '-', color = 'black', linewidth = 1.5, label = 'Labrador basin')
plot(time_rcp, salt_lab_rcp, ':', color = 'gray', linewidth = 1.5)
graph_irm	= plot(time_control, salt_irm_control, '-', color = 'red', linewidth = 1.5, label = 'Irminger basin')
plot(time_rcp, salt_irm_rcp, ':', color = 'firebrick', linewidth = 1.5)
graph_ice	= plot(time_control, salt_ice_control, '-', color = 'blue', linewidth = 1.5, label = 'Iceland basin')
plot(time_rcp, salt_ice_rcp, ':', color = 'cyan', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Salinity (g kg$^{-1}$)')
ax.set_ylim(34.7, 35.3)
ax.set_xlim(1600, 2100)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_lab + graph_irm + graph_ice

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc= (0.78, 0.9), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) Salinity NA regions, LR-CESM')

ax3 	= fig.add_axes([0.555, 0.13, 0.3, 0.3], projection = ccrs.Stereographic(central_longitude=-40, central_latitude=40))

CS      = ax3.contourf(lon_ocn, lat_ocn, salt_1000_3000_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), fraction=0.034, pad=0.04)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)', fontsize = 6.5)

gl = ax3.gridlines(draw_labels=False)
ax3.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax3.coastlines('50m')
ax3.add_feature(cfeature.LAND, zorder=0)
ax3.set_title('Salinity (1000 - 3000 m) trend, \n PI control (21 - 100)', fontsize = 10)

#Labrador basin
lon_lab	= [-44, -50, -63, -63, -55, -44]
lat_lab	= [60, 65, 65, 53, 53, 60]

ax3.plot(lon_lab, lat_lab, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Irminger basin
lon_irm	= [-44, -29.4, -20, -44, -44]
lat_irm = [60, 60, 66, 66, 60]

ax3.plot(lon_irm, lat_irm, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Iceland basin
lon_ice	= [-29.4, -19, -7, -20, -29.4]
lat_ice = [60, 58, 62, 66, 60]

ax3.plot(lon_ice, lat_ice, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.16, 0.13, 0.3, 0.3], projection = ccrs.Stereographic(central_longitude=-40, central_latitude=40))

CS      = ax2.contourf(lon_ocn, lat_ocn, salt_0_100_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), fraction=0.034, pad=0.04)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)', fontsize = 6.5)

gl = ax2.gridlines(draw_labels=False)
ax2.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('Salinity (0 - 100 m) trend, \n PI control (21 - 100)', fontsize = 10)

#-----------------------------------------------------------------------------------------

ax3 	= fig.add_axes([0.555, 0.13, 0.3, 0.3], projection = ccrs.Stereographic(central_longitude=-40, central_latitude=40))

CS      = ax3.contourf(lon_ocn, lat_ocn, salt_1000_3000_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), fraction=0.034, pad=0.04)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)', fontsize = 6.5)

gl = ax3.gridlines(draw_labels=False)
ax3.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax3.coastlines('50m')
ax3.add_feature(cfeature.LAND, zorder=0)
ax3.set_title('Salinity (1000 - 3000 m) trend, \n PI control (21 - 100)', fontsize = 10)

#Labrador basin
lon_lab	= [-44, -50, -63, -63, -55, -44]
lat_lab	= [60, 65, 65, 53, 53, 60]

ax3.plot(lon_lab, lat_lab, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Irminger basin
lon_irm	= [-44, -29.4, -20, -44, -44]
lat_irm = [60, 60, 66, 66, 60]

ax3.plot(lon_irm, lat_irm, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Iceland basin
lon_ice	= [-29.4, -19, -7, -20, -29.4]
lat_ice = [60, 58, 62, 66, 60]

ax3.plot(lon_ice, lat_ice, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

show()

