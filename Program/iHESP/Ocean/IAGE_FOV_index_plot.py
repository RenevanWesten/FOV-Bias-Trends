#Program plots the ideal age at 34S

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	    = '../../../Data/iHESP/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	age		    = fh.variables['IAGE'][:]	    #Ideal age (yr)
	age_ASW		= fh.variables['IAGE_ASW'][:]	#Ideal age (yr)
	age_AIW		= fh.variables['IAGE_AIW'][:]	#Ideal age (yr)
	age_NADW	= fh.variables['IAGE_NADW'][:]	#Ideal age (yr)
	age_ABW		= fh.variables['IAGE_ABW'][:]	#Ideal age (yr)

	fh.close()

	return time, age, age_ASW, age_AIW, age_NADW, age_ABW

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	
section_name	= 'section_34S'

#Read in the time series of the PI control
time_control_HR, age_control_HR, age_ASW_control_HR, age_AIW_control_HR, age_NADW_control_HR, age_ABW_control_HR	= ReadinData(directory+'Ocean/IAGE_'+section_name+'_PI_control_HR.nc')
time_control_LR, age_control_LR, age_ASW_control_LR, age_AIW_control_LR, age_NADW_control_LR, age_ABW_control_LR	= ReadinData(directory+'Ocean/IAGE_'+section_name+'_PI_control_LR.nc')
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control		= plot(time_control_HR, age_ASW_control_HR, '-r', linewidth = 2.0, label = 'HR-CESM PI control')
graph_control_LR	= plot(time_control_LR, age_ASW_control_LR, '-b', linewidth = 2.0, label = 'LR-CESM PI control')
plot([1, 500], [1, 500], '--k')

ax.set_xlabel('Model year')
ax.set_ylabel('Age (year)')
ax.set_xlim(1, 500)
ax.set_ylim(1, 500)
ax.grid()
ax.set_xticks([1, 100, 200, 300, 400, 500])
ax.set_xticklabels(['1', '100', '200', '300', '400',  '500'])
ax.set_yticks([1, 100, 200, 300, 400, 500])
ax.set_yticklabels(['1', '100', '200', '300', '400',  '500'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_control_LR

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('e) Atlantic Surface Water (ASW)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control		= plot(time_control_HR, age_AIW_control_HR, '-r', linewidth = 2.0, label = 'HR-CESM PI control')
graph_control_LR	= plot(time_control_LR, age_AIW_control_LR, '-b', linewidth = 2.0, label = 'LR-CESM PI control')
plot([1, 500], [1, 500], '--k')

ax.set_xlabel('Model year')
ax.set_ylabel('Age (year)')
ax.set_xlim(1, 500)
ax.set_ylim(1, 500)
ax.grid()
ax.set_xticks([1, 100, 200, 300, 400, 500])
ax.set_xticklabels(['1', '100', '200', '300', '400',  '500'])
ax.set_yticks([1, 100, 200, 300, 400, 500])
ax.set_yticklabels(['1', '100', '200', '300', '400',  '500'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_control_LR

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('f) Antarctic Intermediate Water (AIW)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control		= plot(time_control_HR, age_NADW_control_HR, '-r', linewidth = 2.0, label = 'HR-CESM PI control')
graph_control_LR	= plot(time_control_LR, age_NADW_control_LR, '-b', linewidth = 2.0, label = 'LR-CESM PI control')
plot([1, 500], [1, 500], '--k')

ax.set_xlabel('Model year')
ax.set_ylabel('Age (year)')
ax.set_xlim(1, 500)
ax.set_ylim(1, 500)
ax.grid()
ax.set_xticks([1, 100, 200, 300, 400, 500])
ax.set_xticklabels(['1', '100', '200', '300', '400',  '500'])
ax.set_yticks([1, 100, 200, 300, 400, 500])
ax.set_yticklabels(['1', '100', '200', '300', '400',  '500'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_control_LR

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('g) North Atlantic Deep Water (NADW)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control		= plot(time_control_HR, age_ABW_control_HR, '-r', linewidth = 2.0, label = 'HR-CESM PI control')
graph_control_LR	= plot(time_control_LR, age_ABW_control_LR, '-b', linewidth = 2.0, label = 'LR-CESM PI control')
plot([1, 500], [1, 500], '--k')

ax.set_xlabel('Model year')
ax.set_ylabel('Age (year)')
ax.set_xlim(1, 500)
ax.set_ylim(1, 500)
ax.grid()
ax.set_xticks([1, 100, 200, 300, 400, 500])
ax.set_xticklabels(['1', '100', '200', '300', '400',  '500'])
ax.set_yticks([1, 100, 200, 300, 400, 500])
ax.set_yticklabels(['1', '100', '200', '300', '400',  '500'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control + graph_control_LR

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('h) Antarctic Bottom Water (ABW)')

show()


