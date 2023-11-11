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

ax.set_title('d) Salinity NA regions (1000 - 3000 m), LR-CESM')

show()

