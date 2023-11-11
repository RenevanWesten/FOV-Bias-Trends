#Program plots the FOV and FAZ at 34S

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
	FOV		    = fh.variables['F_OV'][:]	#Freshwater transport (Sv)

	fh.close()

	return time, FOV

def ReadinDataGyre(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV_gyre	= fh.variables['F_gyre'][:]	#Freshwater transport (Sv)

	fh.close()

	return time, FOV_gyre

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

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, FOV_control	     = ReadinData(directory+'Ocean/FOV_index_section_34S_PI_control_LR.nc')
time_control1, FW_gyre_control	 = ReadinDataGyre(directory+'Ocean/FW_gyre_section_34S_PI_control_LR.nc')
time_rcp, FOV_rcp		         = ReadinData(directory+'Ocean/FOV_index_section_34S_RCP_LR.nc')
time_rcp, FW_gyre_rcp		     = ReadinDataGyre(directory+'Ocean/FW_gyre_section_34S_RCP_LR.nc')
time_rean, FOV_rean		         = ReadinData(directory_reanalysis+'Ocean/FOV_index_section_34S.nc')
time_rean, FW_gyre_rean	     = ReadinDataGyre(directory_reanalysis+'Ocean/FW_gyre_section_34S.nc')
time_control			         += 1599
#-----------------------------------------------------------------------------------------

print('Trend model year:', int(time_rcp[150]),'-', int(time_rcp[250]))
print()
trend, error, p_value = SignificantTrend(time_rcp[150:251], FOV_rcp[150:251])
print('FovS:', trend*100, '('+str(p_value)+')')
print()
trend, error, p_value = SignificantTrend(time_rcp[150:251], FW_gyre_rcp[150:251])
print('FazS:', trend*100, '('+str(p_value)+')')
print()

#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_control	= ax.plot(time_control, FOV_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	    = ax.plot(time_rcp, FOV_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	    = ax.plot(time_rean, FOV_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.35, 0.35)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) $F_{\mathrm{ovS}}$, LR-CESM')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control  = ax.plot(time_control, FW_gyre_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	   = ax.plot(time_rcp, FW_gyre_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	   = ax.plot(time_rean, FW_gyre_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.1, 0.6)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('d) $F_{\mathrm{azS}}$, LR-CESM')

show()


