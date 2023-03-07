#Program plots the salinity in the Indian Ocean

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
directory	   = '../../../Data/iHESP/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	salt		= fh.variables['SALT_IO'][:]	#Salinity indian Ocean (g/kg)

	fh.close()

	return time, salt

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
depth_max	= 100

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, salt_control	= ReadinData(directory+'Ocean/SALT_Indian_Ocean_'+str(depth_min)+'_'+str(depth_max)+'m_PI_control_HR.nc')
time_rcp, salt_rcp		    = ReadinData(directory+'Ocean/SALT_Indian_Ocean_'+str(depth_min)+'_'+str(depth_max)+'m_RCP_HR.nc')
time_control			    += 1599

#-----------------------------------------------------------------------------------------

print('Trend model year:', int(time_rcp[150]),'-', int(time_rcp[250]))
print()
trend, error, p_value = SignificantTrend(time_rcp[150:251], salt_rcp[150:251])
print('SALT:', trend*100, '('+str(p_value)+')')
print()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory+'Ocean/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_21-100_PI_control_HR.nc', 'r')

lon_ocn		= fh.variables['lon'][:] 	
lat_ocn		= fh.variables['lat'][:] 		
salt_trend	= fh.variables['SALT_trend'][:] 	

fh.close()

#-----------------------------------------------------------------------------------------

fh 		= netcdf.Dataset(directory+'Data/PI_control_HR/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m/iHESP_data_year_0021.nc', 'r')

salt_year_21	= fh.variables['SALT'][:]

fh.close()

fh 		= netcdf.Dataset(directory+'Data/PI_control_HR/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m/iHESP_data_year_0100.nc', 'r')

salt_year_100	= fh.variables['SALT'][:]

fh.close()

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control	= ax.plot(time_control, salt_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= ax.plot(time_rcp, salt_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')

ax.set_xlabel('Model year')
ax.set_ylabel('Salinity (g kg$^{-1}$)')
ax.set_xlim(1600, 2100)
ax.set_ylim(34, 36)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])
ax.set_yticks(np.arange(34, 36.01, 0.5))
ax.set_yticklabels([34.0, 34.5, 35.0, 35.5, 36.0])


graphs	      = graph_control + graph_rcp

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Indian Ocean (0 - 100 m) salinity, HR-CESM')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.20, 0.45, 0.45, 0.30], projection = ccrs.PlateCarree())

CS      = ax2.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())
ax2.contour(lon_ocn, lat_ocn, salt_year_21, levels = [35], linewidths = 0.5, colors = 'k', transform=ccrs.PlateCarree())
ax2.contour(lon_ocn, lat_ocn, salt_year_100, levels = [35], linewidths = 0.5, colors = 'r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.01, 1), fraction=0.021, pad=0.04)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)', fontsize = 6.5)

gl = ax2.gridlines(draw_labels=False)
ax2.set_extent([-80, 130, -70, 25], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('Salinity (0 - 100 m) trend, PI control (21 - 100)', fontsize = 10)


show()


