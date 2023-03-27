#Program plots the FOV at 34S and 60N

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	            = '../../../Data/iHESP/'
directory_reanalysis	= '../../../Data/Reanalysis/'
def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time	= fh.variables['time'][:]		
	FOV		= fh.variables['F_OV'][:]	#Fresh water

	fh.close()

	return time, FOV

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, FOV_34S_control	  = ReadinData(directory+'Ocean/FOV_index_section_34S_PI_control_LR.nc')
time_control, FOV_60N_control	  = ReadinData(directory+'Ocean/FOV_index_section_60N_PI_control_LR.nc')
time_rcp, FOV_34S_rcp		      = ReadinData(directory+'Ocean/FOV_index_section_34S_RCP_LR.nc')
time_rcp, FOV_60N_rcp		      = ReadinData(directory+'Ocean/FOV_index_section_60N_RCP_LR.nc')
time_control			          += 1599
#-----------------------------------------------------------------------------------------

time_rean, FOV_34S_rean = ReadinData(directory_reanalysis+'Ocean/FOV_index_section_34S.nc')
time_rean, FOV_60N_rean = ReadinData(directory_reanalysis+'Ocean/FOV_index_section_60N.nc')

print('FOV RCP LR-CESM (60N): ', np.mean(FOV_60N_rcp[144:171]))
print('FOV Reanalysis (60N): ', np.mean(FOV_60N_rean))
     
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control_34S	= ax.plot(time_control, FOV_34S_control, '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$')
graph_control_60N	= ax.plot(time_control, FOV_60N_control, '-b', linewidth = 1.5, label = '$F_{\mathrm{ovN}}$')
graph_control_conver	= ax.plot(time_control, FOV_34S_control - FOV_60N_control, '-r', linewidth = 1.5, label = '$\Delta F_{\mathrm{ov}}$')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.3, 0.3)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300', '400',  '500'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      = graph_control_34S + graph_control_60N + graph_control_conver

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('b) Freshwater convergence, LR-CESM, PI control')
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_rcp_34S		= ax.plot(time_rcp, FOV_34S_rcp, '-k', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$')
graph_rcp_60N		= ax.plot(time_rcp, FOV_60N_rcp, '-b', linewidth = 1.5, label = '$F_{\mathrm{ovN}}$')
graph_rcp_conver	= ax.plot(time_rcp, FOV_34S_rcp - FOV_60N_rcp, '-r', linewidth = 1.5, label = '$\Delta F_{\mathrm{ov}}$')

graph_rean_34S		= ax.plot(time_rean, FOV_34S_rean, '-', color = 'gray', linewidth = 1.5, label = '$F_{\mathrm{ovS}}$, Reanalysis')
graph_rean_60N		= ax.plot(time_rean, FOV_60N_rean, '-', color = 'cyan', linewidth = 1.5, label = '$F_{\mathrm{ovN}}$, Reanalysis')
graph_rean_conver	= ax.plot(time_rean, FOV_34S_rean - FOV_60N_rean, '-', color = 'firebrick', linewidth = 1.5, label = '$\Delta F_{\mathrm{ov}}$, Reanalysis')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.3, 0.3)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['', '', '', '1900', '2000',  '2100'])

ax.fill_between([1600, 2100], -0.28, -0.05, alpha=0.25, edgecolor='orange', facecolor='orange')

graphs	      	= graph_rcp_34S + graph_rcp_60N + graph_rcp_conver
legend_labels 	= [l.get_label() for l in graphs]
legend_1	= ax.legend(graphs, legend_labels, loc='upper left', ncol=1, framealpha = 1.0, numpoints = 1)


graphs	      	= graph_rean_34S + graph_rean_60N + graph_rean_conver
legend_labels 	= [l.get_label() for l in graphs]
legend_2	= ax.legend(graphs, legend_labels, loc = 'lower left', ncol=1, framealpha = 1.0, numpoints = 1)
ax.add_artist(legend_1)

ax.set_title('d) Freshwater convergence, LR-CESM, Hist/RCP8.5')

show()
