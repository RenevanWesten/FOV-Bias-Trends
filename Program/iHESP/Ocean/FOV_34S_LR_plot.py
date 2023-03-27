#Program plots the FOV at 34S

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
directory	= '../../../Data/iHESP/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Freshwater transport (Sv)

	fh.close()

	return time, FOV

def ReadinDataGyre(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV_gyre	= fh.variables['F_gyre'][:]	#Freshwater transport (Sv)

	fh.close()

	return time, FOV_gyre

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, FOV_control	     = ReadinData(directory+'Ocean/FOV_index_section_34S_PI_control_LR.nc')
time_control1, FOV_gyre_control	 = ReadinDataGyre(directory+'Ocean/FOV_gyre_section_34S_PI_control_LR.nc')
time_rcp, FOV_rcp		         = ReadinData(directory+'Ocean/FOV_index_section_34S_RCP_LR.nc')
time_rcp, FOV_gyre_rcp		     = ReadinDataGyre(directory+'Ocean/FOV_gyre_section_34S_RCP_LR.nc')
time_control			         += 1599
#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Atmosphere/PREC_EVAP_SO_trend_year_21-100_PI_control_LR.nc', 'r')

lon_atm		= fh.variables['lon'][:] 	
lat_atm		= fh.variables['lat'][:] 		
P_E_trend	= fh.variables['P_E_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control	= ax.plot(time_control, FOV_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= ax.plot(time_rcp, FOV_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.3, 0.3)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_rcp

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) $F_{\mathrm{ovS}}$, LR-CESM')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.30, 0.13, 0.45, 0.30], projection = ccrs.PlateCarree())

CS      = ax2.contourf(lon_atm, lat_atm, P_E_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.1, 1), fraction=0.021, pad=0.04)
cbar.set_label('P-E trend (mm day$^{-1}$ per century)', fontsize = 6.5)

gl = ax2.gridlines(draw_labels=False)
ax2.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('P-E trend, PI control (21 - 100)', fontsize = 10)

x_1	= np.arange(-60, 20.01, 0.1)
y_1	= np.zeros(len(x_1)) - 34
y_2	= np.arange(-37, -30.99, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(-37, -30.99, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

ax3 	= fig.add_axes([0.23, 0.65, 0.25, 0.18])
ax3.plot(time_control, FOV_gyre_control, '-k', linewidth = 1.0)
ax3.plot(time_rcp, FOV_gyre_rcp, '-r', linewidth = 1.0)

ax3.set_ylabel('Freshwater transport (Sv)', fontsize = 6.5)
ax3.set_xlim(1600, 2100)
ax3.set_ylim(0, 0.6)
ax3.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax3.set_xticklabels([])
ax3.set_title('Azonal (gyre) component', fontsize = 10)
ax3.grid()

show()


