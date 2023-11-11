#Program plots the AMOC strength at 26N

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
from scipy import stats

#Making pathway to folder with all data
directory	            = '../../../Data/iHESP/'
directory_reanalysis	= '../../../Data/Reanalysis/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#AMOC strength (Sv)

	fh.close()

	return time, transport

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric salterature time series
	time		= np.arange(1, len(time) + 1)

	#Determine the detrended time series
	trend, base 	= polyfit(time, data, 1)
	data_res	= data - ((trend * time) + base)

	#Effective sample size, based on the lag-1 correlation
	corr_1		= np.corrcoef(data_res[:-1], data_res[1:])[0, 1]
	N_eff		= int(len(time) * (1.0 - corr_1) / (1.0 + corr_1))

	#Determine the variance of the anomalies
	data_var	= np.sum(data_res**2.0) / (N_eff - 2.0)

	#Determine the standard error
	standard_error	=  np.sqrt(data_var) / np.sqrt(np.sum((time - np.mean(time))**2.0))

	#Determine the Student-T value
	t_value		= trend / standard_error

	#Get the significance levels and the corresponding critical values (two-sided)
	sig_levels 	= np.arange(50, 100, 0.5) / 100.0
	t_crit 		= stats.t.ppf((1.0 + sig_levels) / 2.0, N_eff - 2)

	#Get the indices where the significance is exceeding the critical values
	sig_index	= np.where(fabs(t_value) > t_crit)[0]
	significant	= 0.0

	if len(sig_index) > 0:
		#If there are significance values, take the highest significant level
		significant = sig_levels[sig_index[-1]]

	return trend, np.sqrt(standard_error), significant


#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	
depth_min	= 0
depth_max	= 1000

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, transport_control	= ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_PI_control_HR.nc')
time_rcp, transport_rcp		    = ReadinData(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP_HR.nc')
time_rean, transport_rean       = ReadinData(directory_reanalysis+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')
time_control			       += 1599
#-----------------------------------------------------------------------------------------

print('Trend model year:', int(time_control[109]-1599),'-', int(time_control[-2]-1599))
print()
trend, error, p_value = SignificantTrend(time_control[109:-1], transport_control[109:-1])
print('AMOC:', trend*100, '('+str(p_value)+')')
print()

print('Trend model year:', int(time_rcp[150]),'-', int(time_rcp[250]))
print()
trend, error, p_value = SignificantTrend(time_rcp[150:251], transport_rcp[150:251])
print('AMOC:', trend*100, '('+str(p_value)+')')
print()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/SST_NA_trend_year_21-100_PI_control_HR.nc', 'r')

lon_ocn		= fh.variables['TLONG'][:] 	
lat_ocn		= fh.variables['TLAT'][:] 		
temp_trend	= fh.variables['TEMP_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([1600, 2100], 16, 19, alpha=0.25, edgecolor='orange', facecolor='orange')

graph_control	= plot(time_control, transport_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= plot(time_rcp, transport_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	= plot(time_rean, transport_rean, '-c', linewidth = 1.5, label = 'Reanalysis')


ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(0, 25)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) AMOC strength at 26$^{\circ}$N, HR-CESM')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.15, 0.145, 0.35, 0.35], projection = ccrs.Stereographic(central_longitude=-40, central_latitude=40))

CS      = ax2.contourf(lon_ocn, lat_ocn, temp_trend, levels = np.arange(-2, 2.1, 0.1), extend = 'both', cmap = 'RdBu_r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-2, 2.1, 1), fraction=0.034, pad=0.04)
cbar.set_label('SST trend ($^{\circ}$C per century)', fontsize = 10)

gl = ax2.gridlines(draw_labels=False)
ax2.set_extent([-81, 1, 23, 78], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('SST trend, PI control (21 - 100)', fontsize = 10)

x_1	= np.arange(-81, -9.99, 0.1)
y_1	= np.zeros(len(x_1)) + 26.0
y_2	= np.arange(24, 28.01, 0.1)
x_2	= np.zeros(len(y_2)) + x_1[0]
y_3	= np.arange(24, 28.01, 0.1)
x_3	= np.zeros(len(y_3)) + x_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

show()
