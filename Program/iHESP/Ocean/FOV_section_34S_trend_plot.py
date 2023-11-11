#Program plots the salinity trends at 34S

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
directory	   = '../../../Data/iHESP/'

def ReadinData(filename, lat_index, depth_min_index, depth_max_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon_u 		= fh.variables['ULONG'][:]							                #Longitude
	lat_u 		= fh.variables['ULAT'][lat_index]							        #Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 			    #Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 				#Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][lat_index] 		                            #Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index] #Meridional velocity (m/s)

	lon_t 		= fh.variables['TLONG'][:]							                #Longitude
	lat_t 		= fh.variables['TLAT'][lat_index:lat_index+2]						#Latitude 
	grid_x_t	= fh.variables['DXT'][lat_index:lat_index+2] 		                   		#Zonal grid cell length (m)
	salt		= fh.variables['SALT'][depth_min_index:depth_max_index, lat_index:lat_index+2] 	#Salinity (g / kg)

	fh.close()
    
	return lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt
   
def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric temperature time series
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

depth_min 	= 0
depth_max	= 6000

experiment  = 'RCP_HR'
#experiment  = 'RCP_LR'

year_start	= 2000
year_end	= 2100


#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/'+experiment+'/FOV_section_34S/iHESP_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-7:-3]	
	year  = int(date[0:4])

	time[year_i] = year

time_start	= (np.abs(time - year_start)).argmin()
time_end	= (np.abs(time - year_end)).argmin() + 1
time		= time[time_start:time_end]
files		= files[time_start:time_end]

year_start	= int(time[0])
year_end	= int(time[-1])

#-----------------------------------------------------------------------------------------


#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat_u 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
depth_t 	= fh.variables['HT'][:] 	#Depth at t-grid (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lat_index		    = 1
depth_u		        = depth_u[lat_index]
depth_t		        = depth_t[lat_index:lat_index+2]

#-----------------------------------------------------------------------------------------

#Determine the section length per depth layer
lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt 	= ReadinData(files[0], lat_index, depth_min_index, depth_max_index)

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt 	= ReadinData(files[0], lat_index, depth_min_index, depth_max_index)
layer_field_u	= ma.masked_all((len(depth), len(lon_u)))
layer_field_t	= ma.masked_all((len(depth), 2, len(lon_t)))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = v_vel[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff_u	= np.sum(layer_field_u, axis = 0) - depth_u

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff_u	= ma.masked_where(depth_diff_u < 0, depth_diff_u)
	depth_diff_u	= depth_diff_u.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field_u[depth_i]	= layer_field_u[depth_i] - depth_diff_u

#Normalise layer field per layer
layer_field_u_norm  = ma.masked_all(shape(layer_field_u))
grid_x_u_norm	    = ma.masked_all((len(depth), len(lon_u)))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_norm[depth_i]		= layer_field_u[depth_i] / np.sum(layer_field_u[depth_i])
    
	#Normalise the length
	grid_x_u_depth          		= ma.masked_array(grid_x_u, mask = v_vel[depth_i].mask)
	grid_x_u_norm[depth_i]  		= grid_x_u_depth / np.sum(grid_x_u_depth)

#-----------------------------------------------------------------------------------------

v_vel_year	= ma.masked_all((len(files), len(depth), len(lon_u)))
vel_mer_year= ma.masked_all((len(files), len(depth)))
salt_year	= ma.masked_all((len(files), len(depth), len(lon_t)))

for year_i in range(len(files)):

	print(year_i)
	lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt 	= ReadinData(files[year_i], lat_index, depth_min_index, depth_max_index)

	#Determine the time mean over the months of choice
	v_vel_year[year_i]	= v_vel * 100.0
	salt_year[year_i]	= salt[:, 0]

	#Determine the meridional transport
	transport	= v_vel * layer_field_u * grid_x_u

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_x_u)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = v_vel - vel_barotropic

	#Save the meridional baroclinic transport (cm/s)
	vel_mer_year[year_i]	= np.sum(vel_baroclinic * grid_x_u_norm, axis = 1) * 100.0


#-----------------------------------------------------------------------------------------

v_vel_trend	    = ma.masked_all((len(depth), len(lon_u)))
v_vel_trend_sig	= ma.masked_all((len(depth), len(lon_u)))
vel_mer_trend	= ma.masked_all(len(depth))
salt_trend	    = ma.masked_all((len(depth), len(lon_t)))
salt_trend_sig	= ma.masked_all((len(depth), len(lon_t)))

for depth_i in range(len(depth)):
	print(depth_i)
	for lon_i in range(len(lon_u)):

		if v_vel_year[0, depth_i, lon_i] is ma.masked:
			#If masked elements, skip
			continue

		if np.all(v_vel_year[:, depth_i, lon_i] == 0.0):
			#All zeros
			continue
	
		trend, error, sig	= SignificantTrend(np.arange(len(files)), v_vel_year[:, depth_i, lon_i])

		#Determine trend per century
		v_vel_trend[depth_i, lon_i]	= trend * 100.0

		#Save the significance
		v_vel_trend_sig[depth_i, lon_i]	= sig

	for lon_i in range(len(lon_t)):

		if salt_year[0, depth_i, lon_i] is ma.masked:
			#If masked elements, skip
			continue

		if np.all(salt_year[:, depth_i, lon_i] == 0.0):
			#All zeros
			continue
	
		trend, error, sig	= SignificantTrend(np.arange(len(files)), salt_year[:, depth_i, lon_i])

		#Determine trend per century
		salt_trend[depth_i, lon_i]	= trend * 100.0

		#Save the significance
		salt_trend_sig[depth_i, lon_i]	= sig

	if vel_mer_year[0, depth_i] is ma.masked:
		continue

	trend, error, sig	= SignificantTrend(np.arange(len(files)), vel_mer_year[:, depth_i])

	#Determine trend per century
	vel_mer_trend[depth_i]	= trend * 100.0

#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

v_vel_trend_sig		= ma.masked_where(v_vel_trend_sig < 0.95, v_vel_trend_sig)
salt_trend_sig		= ma.masked_where(salt_trend_sig < 0.95, salt_trend_sig)
#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-60, 20], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lon_t, depth, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(-1, 1.01, 1))
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)')

ax.set_xlim(-60, 20)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')

CS = contourf(lon_t, depth, salt_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '\\')
CS = contourf(lon_t, depth, salt_trend_sig, levels = np.arange(-1, 1.01, 0.05), extend = 'both', colors='none', linewidth = 2.0, hatches= '//')	

ax.set_xticklabels(['60$^{\circ}$W', '50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

#-----------------------------------------------------------------------------------------

cNorm  		= colors.Normalize(vmin=-1, vmax= 1) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='RdBu_r') 	#Using colormap
color_south 	= scalarMap.to_rgba(-0.5)
color_north 	= scalarMap.to_rgba(0.5)

print(vel_mer_trend)

ax2 = fig.add_axes([0.125, 0.11, 0.05, 0.768])

ax2.plot(vel_mer_trend, depth, '-k', linewidth = 2.0)
ax2.set_xlim(-0.2, 0.2)
ax2.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)

ax2.axvline(x = 0, linestyle = '--', color = 'k')
ax2.fill_betweenx(depth, vel_mer_trend, where = vel_mer_trend >= 0.0, color = color_north, alpha = 0.50)	
ax2.fill_betweenx(depth, vel_mer_trend, where = vel_mer_trend <= 0.0, color = color_south, alpha = 0.50)

labels =  ax.get_yticks()

for label_i in labels:
	ax2.axhline(y = label_i, color = 'gray', linestyle = ':', alpha = 0.5)
	
ax2.set_xticklabels([])
ax2.set_yticklabels([])

if experiment  	== 'RCP_HR':
	ax.set_title('e) Salinity trend at 34$^{\circ}$S, HR-CESM, Hist/RCP8.5 (2000 - 2100)')

if experiment  	== 'RCP_LR':
	ax.set_title('f) Salinity trend at 34$^{\circ}$S, LR-CESM, Hist/RCP8.5 (2000 - 2100)')

show()
#-----------------------------------------------------------------------------------------

