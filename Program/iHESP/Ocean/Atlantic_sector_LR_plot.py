#Program plot the profiles along 50W - 20E

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/iHESP/'


def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	depth		= fh.variables['z_t'][:]	#Depth
	lat_u 		= fh.variables['ULAT'][:]	#Latitude 
	u_vel 		= fh.variables['UVEL'][:] 	#Zonal velocity (m/s)

	lat_t 		= fh.variables['TLAT'][:]	#Latitude 
	temp		= fh.variables['TEMP'][:] 	#Temperature (deg C)
	salt		= fh.variables['SALT'][:] 	#Salinity (g / kg)
	dens		= fh.variables['PD'][:] 	#Potential density (k / m3)

	fh.close()
    
	return lat_u, lat_t, depth, u_vel, temp, salt, dens
			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

experiment  	= 'RCP_LR'

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/'+experiment+'/Atlantic_sector/iHESP_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-7:-3]	
	year  = int(date[0:4])

	time[year_i] = year

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lat_u, lat_t, depth, u_vel, temp, salt, dens = ReadinData(files[0])

#-----------------------------------------------------------------------------------------

#Define empty array's
temp_all		= ma.masked_all((len(time), len(depth), len(lat_t)))
salt_all		= ma.masked_all((len(time), len(depth), len(lat_t)))
u_vel_all		= ma.masked_all((len(time), len(depth), len(lat_u)))
dens_all		= ma.masked_all((len(time), len(depth), len(lat_t)))

for file_i in range(len(files)):
	#Now determine for each month
	print(files[file_i])
	    
	lat_u, lat_t, depth, u_vel, temp, salt, dens = ReadinData(files[file_i])

	#Save the data
	temp_all[file_i]	= temp
	salt_all[file_i]	= salt
	u_vel_all[file_i]	= u_vel
	dens_all[file_i]	= dens
	
#Take the time mean
temp_all	= np.mean(temp_all, axis = 0)
salt_all	= np.mean(salt_all, axis = 0)
u_vel_all	= np.mean(u_vel_all, axis = 0)
dens_all	= np.mean(dens_all, axis = 0)
#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-80, 10], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat_t, depth, salt_all, levels = np.arange(34, 36.01, 0.1), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(34, 36.01, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-71, 1)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-70, 1, 10))
ax.set_xticklabels(['70$^{\circ}$S', '60$^{\circ}$S', '50$^{\circ}$S', '40$^{\circ}$S', '30$^{\circ}$S', '20$^{\circ}$S','10$^{\circ}$S', 'Eq'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_title('c) Salinity, LR-CESM (1994 - 2020)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-80, 10], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat_t, depth, temp_all, levels = np.arange(-2, 20.01, 0.5), extend = 'both', cmap = 'Spectral_r')
cbar	= colorbar(CS, ticks = np.arange(0, 20.01, 5))
cbar.set_label('Temperature ($^{\circ}$C)')

ax.set_xlim(-71, 1)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-70, 1, 10))
ax.set_xticklabels(['70$^{\circ}$S', '60$^{\circ}$S', '50$^{\circ}$S', '40$^{\circ}$S', '30$^{\circ}$S', '20$^{\circ}$S','10$^{\circ}$S', 'Eq'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_title('Temperature, LR-CESM (1994 - 2020)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-80, 10], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat_u, depth, u_vel_all*100, levels = np.arange(-20, 20.01, 1), extend = 'both', cmap = 'RdBu_r')
cbar	= colorbar(CS, ticks = np.arange(-20, 20.01, 5))
cbar.set_label('Zonal velocity (cm s$^{-1}$)')

CS_1	= ax.contour(lat_t, depth, dens_all, levels = [1027], colors = 'k', linewidths = 2)
CS_2	= ax.contour(lat_t, depth, dens_all, levels = [1025, 1025.25, 1025.5, 1025.75, 1026, 1026.25, 1026.5, 1026.75, 1027.25, 1027.5, 1027.75, 1028], colors = 'k', linewidths = 1)
ax.clabel(CS_1, inline=True, fontsize=10, manual = [(-10, 500)])


ax.set_xlim(-71, 1)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-70, 1, 10))
ax.set_xticklabels(['70$^{\circ}$S', '60$^{\circ}$S', '50$^{\circ}$S', '40$^{\circ}$S', '30$^{\circ}$S', '20$^{\circ}$S','10$^{\circ}$S', 'Eq'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_title('f) Zonal velocity and potential density, LR-CESM (1994 - 2020)')

show()
#-----------------------------------------------------------------------------------------



