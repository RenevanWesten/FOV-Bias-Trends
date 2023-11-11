#Program plots the FOV at 34S

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
directory = '../../../Data/POP/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#MOC strength (Sv)
	MOV		    = fh.variables['F_OV'][:]	    #Freshwater
	MOV_ASW		= fh.variables['F_OV_ASW'][:]	#Freshwater
	MOV_AAIW	= fh.variables['F_OV_AAIW'][:]	#Freshwater
	MOV_NADW	= fh.variables['F_OV_NADW'][:]	#Freshwater
	MOV_AABW	= fh.variables['F_OV_AABW'][:]	#Freshwater

	fh.close()

	return time, transport, MOV, MOV_ASW, MOV_AAIW, MOV_NADW, MOV_AABW

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, transport, FOV, FOV_ASW, FOV_AAIW, FOV_NADW, FOV_AABW	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_FOV	= plot(time, FOV, '-k', linewidth = 2.0, label = '$F_{\mathrm{OV}}$')
graph_FOV_ASW	= plot(time, FOV_ASW, '-r', linewidth = 2.0, label = 'ASW')
graph_FOV_AAIW	= plot(time, FOV_AAIW, '-c', linewidth = 2.0, label = 'AAIW')
graph_FOV_NADW	= plot(time, FOV_NADW, '-b', linewidth = 2.0, label = 'NADW')
graph_FOV_AABW	= plot(time, FOV_AABW, '-', color = 'firebrick', linewidth = 2.0, label = 'AABW')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.35, 0.35)
ax.set_xlim(245, 275)
ax.grid()

ax.fill_between([100, 300], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_FOV + graph_FOV_ASW + graph_FOV_AAIW + graph_FOV_NADW + graph_FOV_AABW

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper right', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('f) $F_{\mathrm{ovS}}$, Stand-alone POP')

show()

