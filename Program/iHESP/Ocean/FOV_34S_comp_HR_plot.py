#Program plots the FOV components at 34S

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from scipy import stats

#Making pathway to folder with all data
directory	            = '../../../Data/iHESP/'
directory_reanalysis	= '../../../Data/Reanalysis/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#MOC strength (Sv)
	FOV		    = fh.variables['F_OV'][:]	    #Freshwater
	FOV_ASW		= fh.variables['F_OV_ASW'][:]	#Freshwater
	FOV_AAIW	= fh.variables['F_OV_AAIW'][:]	#Freshwater
	FOV_NADW	= fh.variables['F_OV_NADW'][:]	#Freshwater
	FOV_AABW	= fh.variables['F_OV_AABW'][:]	#Freshwater

	fh.close()

	return time, transport, FOV, FOV_ASW, FOV_AAIW, FOV_NADW, FOV_AABW

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
section_name	= 'section_34S'

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, transport_control, FOV_control, FOV_ASW_control, FOV_AAIW_control, FOV_NADW_control, FOV_AABW_control = ReadinData(directory+'Ocean/FOV_index_'+section_name+'_PI_control_HR.nc')
time_rcp, transport_rcp, FOV_rcp, FOV_ASW_rcp, FOV_AAIW_rcp, FOV_NADW_rcp, FOV_AABW_rcp	= ReadinData(directory+'Ocean/FOV_index_'+section_name+'_RCP_HR.nc')
time_control	 += 1599

#-----------------------------------------------------------------------------------------

time_rean, transport_rean, FOV_rean, FOV_ASW_rean, FOV_AAIW_rean, FOV_NADW_rean, FOV_AABW_rean = ReadinData(directory_reanalysis+'Ocean/FOV_index_'+section_name+'.nc')

#-----------------------------------------------------------------------------------------

print('Trend model year:', int(time_control[229]-1599),'-', int(time_control[-2]-1599))
print()
trend, error, p_value = SignificantTrend(time_control[229:-1], FOV_NADW_control[229:-1])
print('NADW:', trend*100, '('+str(p_value)+')')
print()

#-----------------------------------------------------------------------------------------


trend_1, err, sig	= SignificantTrend(time_rcp[150:], FOV_rcp[150:])
trend_2, err, sig	= SignificantTrend(time_rcp[150:], FOV_ASW_rcp[150:])
trend_3, err, sig	= SignificantTrend(time_rcp[150:], FOV_AAIW_rcp[150:])
trend_4, err, sig	= SignificantTrend(time_rcp[150:], FOV_NADW_rcp[150:])
trend_5, err, sig	= SignificantTrend(time_rcp[150:], FOV_AABW_rcp[150:])

print('Trend model year:', int(time_rcp[150]),'-', int(time_rcp[-1]))
print()
trend, error, p_value = SignificantTrend(time_rcp[150:], FOV_rcp[150:])
print('FOV:', trend*100, '('+str(p_value)+')')
print()

print('ASW trend contribution:', trend_2 / trend_1 * 100, '%')
print('AAIW trend contribution:', trend_3 / trend_1 * 100, '%')
print('NADW trend contribution:', trend_4 / trend_1 * 100, '%')
print('AABW trend contribution:', trend_5 / trend_1 * 100, '%')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

plot(time_control, FOV_control, '-', color = 'black', linewidth = 1.0, alpha = 0.3, label = 'PI control')
plot(time_rcp, FOV_rcp, '-', color = 'red', linewidth = 1.0, alpha = 0.3, label = 'RCP 8.5')
graph_control	= plot(time_control, FOV_ASW_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= plot(time_rcp, FOV_ASW_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	= plot(time_rean, FOV_ASW_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(1600, 2100)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) Atlantic Surface Water (ASW), HR-CESM')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

plot(time_control, FOV_control, '-', color = 'black', linewidth = 1.0, alpha = 0.3, label = 'PI control')
plot(time_rcp, FOV_rcp, '-', color = 'red', linewidth = 1.0, alpha = 0.3, label = 'RCP 8.5')
graph_control	= plot(time_control, FOV_AAIW_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= plot(time_rcp, FOV_AAIW_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	= plot(time_rean, FOV_AAIW_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(1600, 2100)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('c) Antarctic Intermediate Water (AAIW), HR-CESM')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

plot(time_control, FOV_control, '-', color = 'black', linewidth = 1.0, alpha = 0.3, label = 'PI control')
plot(time_rcp, FOV_rcp, '-', color = 'red', linewidth = 1.0, alpha = 0.3, label = 'RCP 8.5')
graph_control	= plot(time_control, FOV_NADW_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= plot(time_rcp, FOV_NADW_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	= plot(time_rean, FOV_NADW_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(1600, 2100)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('e) North Atlantic Deep Water (NADW), HR-CESM')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

plot(time_control, FOV_control, '-', color = 'black', linewidth = 1.0, alpha = 0.3, label = 'PI control')
plot(time_rcp, FOV_rcp, '-', color = 'red', linewidth = 1.0, alpha = 0.3, label = 'RCP 8.5')
graph_control	= plot(time_control, FOV_AABW_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= plot(time_rcp, FOV_AABW_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')
graph_rean	= plot(time_rean, FOV_AABW_rean, '-c', linewidth = 1.5, label = 'Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(1600, 2100)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp + graph_rean

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('g) Antarctic Bottom Water (AABW), HR-CESM')

#-----------------------------------------------------------------------------------------


show()


