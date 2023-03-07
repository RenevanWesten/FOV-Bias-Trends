#Program determines the gyre freshwater transport

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/iHESP/'


def ReadinData(filename, lat_index, depth_min_index, depth_max_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon_u 		= fh.variables['ULONG'][:]							#Longitude
	lat_u 		= fh.variables['ULAT'][lat_index]						#Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 				#Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 				#Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][lat_index] 		                               	#Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index] 		#Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][:]							#Longitude
	lat_t 		= fh.variables['TLAT'][lat_index:lat_index+2]					#Latitude 
	grid_x_t	= fh.variables['DXT'][lat_index:lat_index+2] 		                   	#Zonal grid cell length (m)
	salt		= fh.variables['SALT'][depth_min_index:depth_max_index, lat_index:lat_index+2] 	#Salinity (g / kg)
	

	fh.close()
    
	return lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

experiment  	= 'PI_control_HR'
#experiment  	= 'RCP_HR'
#experiment  	= 'PI_control_LR'
#experiment  	= 'RCP_LR'

section_name	= 'FOV_section_34S'

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/'+experiment+'/'+section_name+'/iHESP_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  = files[year_i][-7:-3]	
	year  = int(date[0:4])
	time[year_i] = year

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat_u 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
depth_t 	= fh.variables['HT'][:] 	#Depth at t-grid (m)
area		= fh.variables['TAREA'][:]	#Area at T-grid
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lat_index		    = 1
depth_u		        = depth_u[lat_index]
depth_t		        = depth_t[lat_index:lat_index+2]
area			    = area[lat_index:lat_index+2]

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt 	= ReadinData(files[0], lat_index, depth_min_index, depth_max_index)
layer_field_u					= ma.masked_all((len(depth), len(lon_u)))
layer_field_t					= ma.masked_all((len(depth), 2, len(lon_t)))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = v_vel[depth_i].mask)
	layer_field_t[depth_i]	= layer[depth_i]
	layer_field_t[depth_i]	= ma.masked_array(layer_field_t[depth_i], mask = salt[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff_u	= np.sum(layer_field_u, axis = 0) - depth_u
	depth_diff_t	= np.sum(layer_field_t, axis = 0) - depth_t

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff_u	= ma.masked_where(depth_diff_u < 0, depth_diff_u)
	depth_diff_u	= depth_diff_u.filled(fill_value = 0.0)
	depth_diff_t	= ma.masked_where(depth_diff_t < 0, depth_diff_t)
	depth_diff_t	= depth_diff_t.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field_u[depth_i]	= layer_field_u[depth_i] - depth_diff_u
	layer_field_t[depth_i]	= layer_field_t[depth_i] - depth_diff_t

#Normalise layer field per layer
volume_t	    = area * layer_field_t
layer_field_u_area  = ma.masked_all((len(depth), len(lon_u)))
layer_field_t_area  = ma.masked_all((len(depth), 2, len(lon_t)))
grid_x_u_norm	    = ma.masked_all((len(depth), len(lon_u)))	
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t)))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_area[depth_i]		= layer_field_u[depth_i] * grid_x_u
	layer_field_t_area[depth_i]		= layer_field_t[depth_i] * grid_x_t
    
	#Normalise the length
	grid_x_u_depth          		= ma.masked_array(layer_field_u[depth_i] * grid_x_u, mask = v_vel[depth_i].mask)
	grid_x_t_depth          		= ma.masked_array(layer_field_t[depth_i] * grid_x_t, mask = salt[depth_i].mask)

	#Normalse per depth layer
	grid_x_u_norm[depth_i]			= grid_x_u_depth / np.sum(grid_x_u_depth)

	for lat_i in range(2):
		#Normalise for each latitude
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_gyre_all		= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt = ReadinData(files[file_i], lat_index, depth_min_index, depth_max_index)

	#Determine the zonal means
	v_vel_zonal     = np.sum(v_vel * grid_x_u_norm, axis = 1)
	salt_zonal      = np.sum(salt * grid_x_t_norm, axis = 2)

	v_vel_prime	= ma.masked_all(np.shape(v_vel))
	salt_prime	= ma.masked_all(np.shape(v_vel))


	for depth_i in range(len(depth)):
		#Determine the differences with respect to the zonal mean
		v_vel_prime[depth_i]	= v_vel[depth_i] - v_vel_zonal[depth_i]

		salt_prime_depth	= ma.masked_all((2, len(lon_t)))
		salt_prime_depth[0]	= salt[depth_i, 0] - salt_zonal[depth_i, 0]
		salt_prime_depth[1]	= salt[depth_i, 1] - salt_zonal[depth_i, 1]

		for lon_i in range(len(lon_u)):
			#Now interpolate the S' to the u-grid by taking the volume-averaged mean
			if v_vel[depth_i, lon_i] is ma.masked or layer_field_u[depth_i, lon_i] <= 0.0:
				continue

			#Take the 4 T-grid cells around the current U-cell
			salt_prime_grid	= salt_prime_depth[:, lon_i:lon_i+2]
			volume_t_norm	= volume_t[depth_i, :, lon_i:lon_i+2] / np.sum(volume_t[depth_i, :, lon_i:lon_i+2])

			#Take the mean over the 4 grids
			salt_prime[depth_i, lon_i]	= np.sum(salt_prime_grid*volume_t_norm)

	#Now determine the azonal component (gyre, in Sv)
	transport_gyre_all[file_i]	= (-1.0 / 35.0) * np.sum(v_vel_prime * salt_prime * layer_field_u_area) / 10**6.0
	#-----------------------------------------------------------------------------------------	


#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control	= ax.plot(time, transport_gyre_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(0, 0.6)
ax.grid()

show()

#-----------------------------------------------------------------------------------------
#Do not overwrite the data, the data in the example directory is not the complete time series
#The time series contains all the availble model years
sys.exit()	

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_gyre_'+section_name[4:]+'_'+experiment+'.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('F_gyre', float, ('time'), zlib=True)

fh.variables['F_gyre'].longname 	= 'Freshwater transport by gyre'

fh.variables['time'].units 		    = 'Year'
fh.variables['F_gyre'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['F_gyre'][:] 		    = transport_gyre_all

fh.close()
