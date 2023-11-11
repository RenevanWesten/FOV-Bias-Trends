#Program determines the FOV index

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
	lon_u 		= fh.variables['ULONG'][:]							             #Longitude
	lat_u 		= fh.variables['ULAT'][lat_index]						         #Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 		     #Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 		     #Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][lat_index] 		                               	#Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index] 		#Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][:]							            #Longitude
	lat_t 		= fh.variables['TLAT'][lat_index:lat_index+2]					#Latitude 
	grid_x_t	= fh.variables['DXT'][lat_index:lat_index+2] 		            #Zonal grid cell length (m)
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
layer_field_u_area  = ma.masked_all((len(depth), len(lon_u)))
layer_field_t_area  = ma.masked_all((len(depth), 2, len(lon_t)))
grid_x_t_norm	    = ma.masked_all((len(depth), 2, len(lon_t)))
lat_weight	    = ma.masked_all((len(depth), 2))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_area[depth_i]		= layer_field_u[depth_i] * grid_x_u
	layer_field_t_area[depth_i]		= layer_field_t[depth_i] * grid_x_t
    
	#Normalise the length
	grid_x_t_depth          		= ma.masked_array(grid_x_t, mask = salt[depth_i].mask)

	#Now get the lat weights
	lat_weight[depth_i]			= np.sum(grid_x_t_depth, axis = 1) / np.sum(grid_x_t_depth)

	for lat_i in range(2):
		grid_x_t_norm[depth_i, lat_i]  = grid_x_t_depth[lat_i] / np.sum(grid_x_t_depth[lat_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))
transport_salt_all	= ma.masked_all(len(time))
transport_salt_ASW_all	= ma.masked_all(len(time))
transport_salt_AAIW_all	= ma.masked_all(len(time))
transport_salt_NADW_all	= ma.masked_all(len(time))
transport_salt_AABW_all	= ma.masked_all(len(time))
salt_ASW_all		= ma.masked_all(len(time))
salt_AAIW_all		= ma.masked_all(len(time))
salt_NADW_all		= ma.masked_all(len(time))
salt_AABW_all		= ma.masked_all(len(time))
vel_ASW_all		= ma.masked_all(len(time))
vel_AAIW_all		= ma.masked_all(len(time))
vel_NADW_all		= ma.masked_all(len(time))
vel_AABW_all		= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon_u, lat_u, lon_t, lat_t, depth, layer, grid_x_u, grid_x_t, v_vel, salt = ReadinData(files[file_i], lat_index, depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport	= v_vel * layer_field_u * grid_x_u

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_x_u)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = v_vel - vel_barotropic

	#Determine the zonal means
	salt_zonal      = np.sum(salt * grid_x_t_norm, axis = 2)  - 35.0
	transport_clin	= np.sum(vel_baroclinic * layer_field_u * grid_x_u, axis = 1)

	#-----------------------------------------------------------------------------------------
	#Get the water properties
	water_prop	= ma.masked_all((len(depth), len(lon_u)))
	water_prop_salt	= ma.masked_all((len(depth), 2, len(lon_u)))

	#North Atlantic Deep Water (NADW) has negative meridional velocities
	depth_index_NADW = np.where((depth >= 700) & (transport_clin <= 0))[0][0]

	#Antarctic bottom water (AABW) is directly below the NADW, get the first index
	depth_index_AABW	= np.where((depth >= 3000) & (transport_clin >= 0))[0]

	if len(depth_index_AABW) == 0:
		#Always assume below 4000m depths AABW water mass (in case of no AABW detection)
		depth_index_AABW	= np.where(depth >= 4000)[0][0]
	else:
		depth_index_AABW	= depth_index_AABW[0]

	for depth_i in range(len(depth)):
			
		if depth_i < depth_index_NADW:
			#Surface water
			water_prop[depth_i]	= 1.0
			water_prop_salt[depth_i]= 1.0

		if depth[depth_i] >= 500 and depth_i < depth_index_NADW and section_name != 'MOV_section_26N':
			#Antarctic Intermediate water
			water_prop[depth_i]	= 2.0
			water_prop_salt[depth_i]= 2.0
		
		if depth_i >= depth_index_NADW and depth_i < depth_index_AABW:
			#North Atlantic Deep Water (NADW)
			water_prop[depth_i]	= 3.0	
			water_prop_salt[depth_i]= 3.0	

		if depth_i >= depth_index_AABW:
			#The AABW is defined below the NADW 
			water_prop[depth_i]	= 4.0
			water_prop_salt[depth_i]= 4.0	

	water_prop	= ma.masked_array(water_prop, mask = v_vel.mask)	
	water_prop_salt	= ma.masked_array(water_prop_salt, mask = salt.mask)	
	#-----------------------------------------------------------------------------------------
	area_u_ASW	= ma.masked_where(water_prop != 1.0, layer_field_u_area)
	area_u_AAIW	= ma.masked_where(water_prop != 2.0, layer_field_u_area)
	area_u_NADW	= ma.masked_where(water_prop != 3.0, layer_field_u_area)
	area_u_AABW	= ma.masked_where(water_prop != 4.0, layer_field_u_area)
	area_u_ASW	= area_u_ASW	/ np.sum(area_u_ASW)
	area_u_AAIW	= area_u_AAIW	/ np.sum(area_u_AAIW)
	area_u_NADW	= area_u_NADW	/ np.sum(area_u_NADW)
	area_u_AABW	= area_u_AABW	/ np.sum(area_u_AABW)

	area_t_ASW	= ma.masked_where(water_prop_salt != 1.0, layer_field_t_area)
	area_t_AAIW	= ma.masked_where(water_prop_salt != 2.0, layer_field_t_area)
	area_t_NADW	= ma.masked_where(water_prop_salt != 3.0, layer_field_t_area)
	area_t_AABW	= ma.masked_where(water_prop_salt != 4.0, layer_field_t_area)
	area_t_ASW	= area_t_ASW	/ np.sum(area_t_ASW)
	area_t_AAIW	= area_t_AAIW	/ np.sum(area_t_AAIW)
	area_t_NADW	= area_t_NADW	/ np.sum(area_t_NADW)
	area_t_AABW	= area_t_AABW	/ np.sum(area_t_AABW)

	#Determine the spatial means
	vel_ASW_all[file_i]	= np.sum(vel_baroclinic * area_u_ASW)
	vel_AAIW_all[file_i]	= np.sum(vel_baroclinic * area_u_AAIW)
	vel_NADW_all[file_i]	= np.sum(vel_baroclinic * area_u_NADW)
	vel_AABW_all[file_i]	= np.sum(vel_baroclinic * area_u_AABW)
	salt_ASW_all[file_i]	= np.sum(salt * area_t_ASW)
	salt_AAIW_all[file_i]	= np.sum(salt * area_t_AAIW)
	salt_NADW_all[file_i]	= np.sum(salt * area_t_NADW)
	salt_AABW_all[file_i]	= np.sum(salt * area_t_AABW)

	#Determine the means over the water masses
	transport_ASW	= np.sum(ma.masked_where(water_prop != 1.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
	transport_AAIW	= np.sum(ma.masked_where(water_prop != 2.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
	transport_NADW	= np.sum(ma.masked_where(water_prop != 3.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)
	transport_AABW	= np.sum(ma.masked_where(water_prop!= 4.0, vel_baroclinic * layer_field_u * grid_x_u), axis = 1)

	#Take the mean over the two latitudes for the salt
	salt_zonal      = np.sum(salt_zonal * lat_weight, axis = 1)

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[file_i]		= np.sum(transport) / 1000000.0
	        
	#Determine the total salinity transport
	transport_salt_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 
	transport_salt_ASW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_ASW * salt_zonal) / 1000000.0 
	transport_salt_AAIW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_AAIW * salt_zonal) / 1000000.0 
	transport_salt_NADW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_NADW * salt_zonal) / 1000000.0 
	transport_salt_AABW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_AABW * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_control	= ax.plot(time, transport_salt_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_ylim(-0.3, 0.3)
ax.grid()

show()

#-----------------------------------------------------------------------------------------
#Do not overwrite the data, the data in the example directory is not the complete time series
#The time series contains all the availble model years
sys.exit()	

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_index_'+section_name[4:]+'_'+experiment+'_TEST.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)
fh.createVariable('F_OV', float, ('time'), zlib=True)
fh.createVariable('F_OV_ASW', float, ('time'), zlib=True)
fh.createVariable('F_OV_AAIW', float, ('time'), zlib=True)
fh.createVariable('F_OV_NADW', float, ('time'), zlib=True)
fh.createVariable('F_OV_AABW', float, ('time'), zlib=True)
fh.createVariable('SALT_ASW', float, ('time'), zlib=True)
fh.createVariable('SALT_AAIW', float, ('time'), zlib=True)
fh.createVariable('SALT_NADW', float, ('time'), zlib=True)
fh.createVariable('SALT_AABW', float, ('time'), zlib=True)
fh.createVariable('VVEL_ASW', float, ('time'), zlib=True)
fh.createVariable('VVEL_AAIW', float, ('time'), zlib=True)
fh.createVariable('VVEL_NADW', float, ('time'), zlib=True)
fh.createVariable('VVEL_AABW', float, ('time'), zlib=True)

fh.variables['Transport'].long_name 	= 'Volume transport'
fh.variables['F_OV'].long_name 		    = 'Fresh water transport'
fh.variables['F_OV_ASW'].long_name 	    = 'Fresh water transport (Atlantic Surface Water)'
fh.variables['F_OV_AAIW'].long_name 	= 'Fresh water transport (Antarctic Intermediate Water)'
fh.variables['F_OV_NADW'].long_name 	= 'Fresh water transport (North Atlantic Deep Water)'
fh.variables['F_OV_AABW'].long_name 	= 'Fresh water transport (Antarctic Bottom Water)'
fh.variables['SALT_ASW'].long_name 	    = 'Salinity (Atlantic Surface Water)'
fh.variables['SALT_AAIW'].long_name 	= 'Salinity (Antarctic Intermediate Water)'
fh.variables['SALT_NADW'].long_name 	= 'Salinity (North Atlantic Deep Water)'
fh.variables['SALT_AABW'].long_name 	= 'Salinity (Antarctic Bottom Water)'
fh.variables['VVEL_ASW'].long_name 	    = 'Meridional velocity (Atlantic Surface Water)'
fh.variables['VVEL_AAIW'].long_name 	= 'Meridional velocity (Antarctic Intermediate Water)'
fh.variables['VVEL_NADW'].long_name 	= 'Meridional velocity (North Atlantic Deep Water)'
fh.variables['VVEL_AABW'].long_name 	= 'Meridional velocity (Antarctic Bottom Water)'

fh.variables['time'].units 		    = 'Year'
fh.variables['Transport'].units 	= 'Sv'
fh.variables['F_OV'].units 		    = 'Sv'
fh.variables['F_OV_ASW'].units 		= 'Sv'
fh.variables['F_OV_AAIW'].units 	= 'Sv'
fh.variables['F_OV_NADW'].units 	= 'Sv'
fh.variables['F_OV_AABW'].units 	= 'Sv'
fh.variables['SALT_ASW'].units 		= 'g/kg'
fh.variables['SALT_AAIW'].units 	= 'g/kg'
fh.variables['SALT_NADW'].units 	= 'g/kg'
fh.variables['SALT_AABW'].units 	= 'g/kg'
fh.variables['VVEL_ASW'].units 		= 'cm/s'
fh.variables['VVEL_AAIW'].units 	= 'cm/s'
fh.variables['VVEL_NADW'].units 	= 'cm/s'
fh.variables['VVEL_AABW'].units 	= 'cm/s'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all
fh.variables['F_OV'][:] 		= transport_salt_all
fh.variables['F_OV_ASW'][:] 		= transport_salt_ASW_all
fh.variables['F_OV_AAIW'][:] 		= transport_salt_AAIW_all
fh.variables['F_OV_NADW'][:] 		= transport_salt_NADW_all
fh.variables['F_OV_AABW'][:] 		= transport_salt_AABW_all
fh.variables['SALT_ASW'][:] 		= salt_ASW_all
fh.variables['SALT_AAIW'][:] 		= salt_AAIW_all
fh.variables['SALT_NADW'][:] 		= salt_NADW_all
fh.variables['SALT_AABW'][:] 		= salt_AABW_all
fh.variables['VVEL_ASW'][:] 		= vel_ASW_all * 100.0
fh.variables['VVEL_AAIW'][:] 		= vel_AAIW_all * 100.0
fh.variables['VVEL_NADW'][:] 		= vel_NADW_all * 100.0
fh.variables['VVEL_AABW'][:] 		= vel_AABW_all * 100.0

fh.close()
