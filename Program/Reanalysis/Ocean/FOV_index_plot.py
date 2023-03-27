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
directory = '../../../Data/Reanalysis/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	transport	= fh.variables['Transport'][:]	#MOC strength (Sv)
	MOV		    = fh.variables['F_OV'][:]	    #Freshwater
	MOV_ASW		= fh.variables['F_OV_ASW'][:]	#Freshwater
	MOV_AIW		= fh.variables['F_OV_AIW'][:]	#Freshwater
	MOV_NADW	= fh.variables['F_OV_NADW'][:]	#Freshwater
	MOV_ABW		= fh.variables['F_OV_ABW'][:]	#Freshwater

	fh.close()

	return time, transport, MOV, MOV_ASW, MOV_AIW, MOV_NADW, MOV_ABW

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

time, transport, FOV, FOV_ASW, FOV_AIW, FOV_NADW, FOV_ABW	= ReadinData(directory+'Ocean/FOV_index_section_34S.nc')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_FOV	= plot(time, FOV, '-k', linewidth = 2.0, label = '$F_{\mathrm{ovS}}$')
graph_FOV_ASW	= plot(time, FOV_ASW, '-r', linewidth = 2.0, label = 'ASW')
graph_FOV_AIW	= plot(time, FOV_AIW, '-c', linewidth = 2.0, label = 'AIW')
graph_FOV_NADW	= plot(time, FOV_NADW, '-b', linewidth = 2.0, label = 'NADW')
graph_FOV_ABW	= plot(time, FOV_ABW, '-', color = 'firebrick', linewidth = 2.0, label = 'ABW')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.3, 0.3)
ax.set_xlim(1990, 2020)
ax.grid()

ax.fill_between([1990, 2021], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_FOV + graph_FOV_ASW + graph_FOV_AIW + graph_FOV_NADW + graph_FOV_ABW

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc=(0.8, 0.75), ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('e) $F_{\mathrm{ovS}}$, Reanalysis')

show()

