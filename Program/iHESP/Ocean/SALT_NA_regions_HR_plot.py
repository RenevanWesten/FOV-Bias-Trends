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
time_control, salt_lab_control, salt_irm_control, salt_ice_control	= ReadinData(directory+'Ocean/SALT_NA_regions_PI_control_HR.nc')
time_rcp, salt_lab_rcp, salt_irm_rcp, salt_ice_rcp			= ReadinData(directory+'Ocean/SALT_NA_regions_RCP_HR.nc')
time_control								+= 1599
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
ax.legend(graphs, legend_labels, loc= 'lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Salinity NA regions (1000 - 3000 m), HR-CESM')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.200, 0.13, 0.45, 0.45], projection = ccrs.NearsidePerspective(-30, 60, 2000000))

ax2.gridlines(zorder=10)
ax2.add_feature(cfeature.LAND, zorder=10)
ax2.coastlines()
ax2.set_global()

#Labrador basin
lon_lab	= [-44, -50, -63, -63, -55, -44]
lat_lab	= [60, 65, 65, 53, 53, 60]

lon_lab_1 	= np.arange(-50, -44+0.01, 0.1)[::-1]
lat_lab_1 	= np.linspace(60, 65, len(lon_lab_1), endpoint = True)
lon_lab_2 	= np.arange(-62.9, -50, 0.1)[::-1]
lat_lab_2	= np.zeros(len(lon_lab_2)) + 65.0
lat_lab_3 	= np.arange(53, 65.1, 0.1)[::-1]
lon_lab_3	= np.zeros(len(lat_lab_3))-63
lon_lab_4 	= np.arange(-62.9, -54.9, 0.1)
lat_lab_4	= np.zeros(len(lon_lab_4)) + 53
lon_lab_5 	= np.arange(-55.1, -44, 0.1)
lat_lab_5 	= np.linspace(53, 60, len(lon_lab_5), endpoint = True)

lon_lab		= np.append(lon_lab_1, lon_lab_2)
lon_lab		= np.append(lon_lab, lon_lab_3)
lon_lab		= np.append(lon_lab, lon_lab_4)
lon_lab		= np.append(lon_lab, lon_lab_5)
lat_lab		= np.append(lat_lab_1, lat_lab_2)
lat_lab		= np.append(lat_lab, lat_lab_3)
lat_lab		= np.append(lat_lab, lat_lab_4)
lat_lab		= np.append(lat_lab, lat_lab_5)

ax2.plot(lon_lab, lat_lab, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Irminger basin
lon_irm	= [-44, -29.4, -20, -44, -44]
lat_irm = [60, 60, 66, 66, 60]

lon_irm_1 	= np.arange(-44, -29.4, 0.1)
lat_irm_1 	= np.zeros(len(lon_irm_1))+60
lon_irm_2 	= np.arange(-29.4, -19.99, 0.1)
lat_irm_2	= np.linspace(60, 66, len(lon_irm_2), endpoint = True)
lon_irm_3 	= np.arange(-44, -20, 0.1)[::-1]
lat_irm_3	= np.zeros(len(lon_irm_3))+66
lat_irm_4	= np.arange(60.1, 66, 0.1)[::-1]
lon_irm_4	= np.zeros(len(lat_irm_4))-44

lon_irm		= np.append(lon_irm_1, lon_irm_2)
lon_irm		= np.append(lon_irm, lon_irm_3)
lon_irm		= np.append(lon_irm, lon_irm_4)
lat_irm		= np.append(lat_irm_1, lat_irm_2)
lat_irm		= np.append(lat_irm, lat_irm_3)
lat_irm		= np.append(lat_irm, lat_irm_4)

ax2.plot(lon_irm, lat_irm, '-r', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#Iceland basin
lon_ice	= [-29.4, -19, -7, -20, -29.4]
lat_ice = [60, 58, 62, 66, 60]

lon_ice_1 	= np.arange(-29.4, -18.9, 0.1)
lat_ice_1 	= np.linspace(60, 58, len(lon_ice_1), endpoint = True)
lon_ice_2 	= np.arange(-19, -6.9, 0.1)
lat_ice_2 	= np.linspace(58, 62, len(lon_ice_2), endpoint = True)
lon_ice_3 	= np.arange(-20, -6.9, 0.1)[::-1]
lat_ice_3 	= np.linspace(62, 66, len(lon_ice_3), endpoint = True)
lon_ice_4 	= np.arange(-29.4, -19.9, 0.1)[::-1]
lat_ice_4	= np.linspace(66, 60, len(lon_ice_4), endpoint = True)

lon_ice		= np.append(lon_ice_1, lon_ice_2)
lon_ice		= np.append(lon_ice, lon_ice_3)
lon_ice		= np.append(lon_ice, lon_ice_4)
lat_ice		= np.append(lat_ice_1, lat_ice_2)
lat_ice		= np.append(lat_ice, lat_ice_3)
lat_ice		= np.append(lat_ice, lat_ice_4)

ax2.plot(lon_ice, lat_ice, '-b', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

show()