#Program determines the FOV index

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf

#Making pathway to folder with all data
directory	= '../../../Data/Reanalysis/'

def ReadinData(filename, lat_index, depth_min_index, depth_max_index, section_name):

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon 		= fh.variables['longitude'][:]						#Longitude
	lat 		= fh.variables['longitude'][lat_index]					#Latitude 
	depth   	= fh.variables['depth'][depth_min_index:depth_max_index] 		#Depth (m)
	v_vel 		= fh.variables['vo'][:, depth_min_index:depth_max_index, lat_index] 	#Meridional velocity (m/s)
	salt		= fh.variables['so'][:, depth_min_index:depth_max_index, lat_index] 	#Salinity (g / kg)
	v_vel_mask 	= fh.variables['vo'][0, depth_min_index:depth_max_index, lat_index] 	#Meridional velocity (m/s)

	fh.close()

	fh		= netcdf.Dataset(directory+'Data/Grid/Grid_'+section_name+'.nc', 'r')

	layer		= fh.variables['layer'][depth_min_index:depth_max_index, lat_index] 	#Layer thickness (m)
	grid_x		= fh.variables['DX'][lat_index] 					#Zonal grid cell length (m)

	fh.close()
    
	#Convert to yearly averaged data
	month_days	= np.asarray([31., 28., 31., 30., 31., 30., 31., 31., 30., 31., 30., 31.])
	month_days	= month_days / np.sum(month_days)

	#Fill the array's with the same dimensions
	month_days_all	= ma.masked_all((len(month_days), len(depth), len(lon)))

	for month_i in range(len(month_days)):
		month_days_all[month_i]		= month_days[month_i]

	#Determine the time mean over the months of choice
	v_vel	= np.sum(v_vel * month_days_all, axis = 0)
	v_vel	= ma.masked_array(v_vel, mask = v_vel_mask.mask)
	salt	= np.sum(salt * month_days_all, axis = 0)
	salt	= ma.masked_array(salt, mask = v_vel_mask.mask)

	return lon, lat, depth, layer, grid_x, v_vel, salt
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

#section_name	= 'FOV_section_34S'
section_name	= 'FOV_section_60N'

#-----------------------------------------------------------------------------------------

files	= glob.glob(directory+'Data/'+section_name+'/Reanalysis_data_year_*.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Define empty array's
time 		= np.zeros(len(files))

for year_i in range(len(files)):
	date  		= files[year_i][-7:-3]	
	year  		= int(date[0:4])
	time[year_i]	= year

#-----------------------------------------------------------------------------------------

if section_name == 'FOV_section_34S':
    lat_FOV = -34
    
if section_name == 'FOV_section_60N':
    lat_FOV = 60

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat 		= fh.variables['latitude'][:]	#Latitude  
depth   	= fh.variables['depth'][:]	#Depth (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lat_index		    = (fabs(lat_FOV - lat)).argmin()

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[0], lat_index, depth_min_index, depth_max_index, section_name)

#Normalise layer field per layer
layer_field_norm  = ma.masked_all(shape(layer_field))
grid_x_norm	  = ma.masked_all((len(depth), len(lon)))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_norm[depth_i]	= layer_field[depth_i] / np.sum(layer_field[depth_i])
    
	#Normalise the length
	grid_x_depth          		= ma.masked_array(grid_x, mask = v_vel[depth_i].mask)
	grid_x_norm[depth_i]  		= grid_x_depth / np.sum(grid_x_depth)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))
transport_salt_all	= ma.masked_all(len(time))
transport_salt_ASW_all	= ma.masked_all(len(time))
transport_salt_AIW_all	= ma.masked_all(len(time))
transport_salt_NADW_all	= ma.masked_all(len(time))
transport_salt_ABW_all	= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon, lat, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[file_i], lat_index, depth_min_index, depth_max_index, section_name)

	#Determine the meridional transport
	transport	= v_vel * layer_field * grid_x

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field * grid_x)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = v_vel - vel_barotropic

	#Determine the zonal means
	salt_zonal      = np.sum(salt * grid_x_norm, axis = 1)  - 35.0
	transport_clin	= np.sum(vel_baroclinic * layer_field * grid_x, axis = 1)

	#-----------------------------------------------------------------------------------------
	#Get the water properties
	water_prop	= ma.masked_all((len(depth), len(lon)))

	#North Atlantic Deep Water (NADW) has negative meridional velocities
	depth_index_NADW = np.where((depth >= 700) & (transport_clin <= 0))[0][0]

	#Antarctic bottom water (ABW) is directly below the NADW, get the first index
	depth_index_ABW	= np.where((depth >= 3000) & (transport_clin >= 0))[0]

	if len(depth_index_ABW) == 0:
		#Assume below 4000m depth the ABW
		depth_index_ABW	= np.where(depth >= 4000)[0][0]
	else:
		depth_index_ABW	= depth_index_ABW[0]

	for depth_i in range(len(depth)):
			
		if depth_i < depth_index_NADW:
			#Surface water
			water_prop[depth_i]	= 1.0

		if depth[depth_i] >= 500 and depth_i < depth_index_NADW and section_name != 'MOV_section_26N':
			#Antarctic Intermediate water
			water_prop[depth_i]	= 2.0
		
		if depth_i >= depth_index_NADW and depth_i < depth_index_ABW:
			#North Atlantic Deep Water (NADW)
			water_prop[depth_i]	= 3.0	

		if depth_i >= depth_index_ABW:
			#The ABW is defined below the NADW 
			water_prop[depth_i]	= 4.0	

	water_prop	= ma.masked_array(water_prop, mask = v_vel.mask)	

	#-----------------------------------------------------------------------------------------

	#Determine the means over the water masses
	transport_ASW	= np.sum(ma.masked_where(water_prop != 1.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_AIW	= np.sum(ma.masked_where(water_prop != 2.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_NADW	= np.sum(ma.masked_where(water_prop != 3.0, vel_baroclinic * layer_field * grid_x), axis = 1)
	transport_ABW	= np.sum(ma.masked_where(water_prop != 4.0, vel_baroclinic * layer_field * grid_x), axis = 1)

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[file_i]		= np.sum(transport) / 1000000.0
	        
	#Determine the total salinity transport
	transport_salt_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_clin * salt_zonal) / 1000000.0 
	transport_salt_ASW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_ASW * salt_zonal) / 1000000.0 
	transport_salt_AIW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_AIW * salt_zonal) / 1000000.0 
	transport_salt_NADW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_NADW * salt_zonal) / 1000000.0 
	transport_salt_ABW_all[file_i]	= (-1.0 / 35.0) * np.sum(transport_ABW * salt_zonal) / 1000000.0 

#-----------------------------------------------------------------------------------------
D
print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FOV_index_'+section_name[4:]+'.nc', 'w')

if section_name == 'FOV_section_60N':
    #Only write the FOV, not the other components
    fh.createDimension('time', len(time))
    
    fh.createVariable('time', float, ('time'), zlib=True)
    fh.createVariable('Transport', float, ('time'), zlib=True)
    fh.createVariable('F_OV', float, ('time'), zlib=True)
    
    fh.variables['Transport'].longname 	= 'Volume transport'
    fh.variables['F_OV'].longname 		= 'Fresh water transport'

    fh.variables['time'].units 		    = 'Year'
    fh.variables['Transport'].units 	= 'Sv'
    fh.variables['F_OV'].units 		    = 'Sv'

    #Writing data to correct variable	
    fh.variables['time'][:]     	  	= time
    fh.variables['Transport'][:]    	= transport_all
    fh.variables['F_OV'][:] 			= transport_salt_all
    
    fh.close()

else:
    #The FOV and its components
    fh.createDimension('time', len(time))
    
    fh.createVariable('time', float, ('time'), zlib=True)
    fh.createVariable('Transport', float, ('time'), zlib=True)
    fh.createVariable('F_OV', float, ('time'), zlib=True)
    fh.createVariable('F_OV_ASW', float, ('time'), zlib=True)
    fh.createVariable('F_OV_AIW', float, ('time'), zlib=True)
    fh.createVariable('F_OV_NADW', float, ('time'), zlib=True)
    fh.createVariable('F_OV_ABW', float, ('time'), zlib=True)
    
    fh.variables['Transport'].longname 	= 'Volume transport'
    fh.variables['F_OV'].longname 		= 'Fresh water transport'
    fh.variables['F_OV_ASW'].longname 	= 'Fresh water transport (Atlantic Surface Water)'
    fh.variables['F_OV_AIW'].longname 	= 'Fresh water transport (Antarctic Intermediate Water)'
    fh.variables['F_OV_NADW'].longname 	= 'Fresh water transport (North Atlantic Deep Water)'
    fh.variables['F_OV_ABW'].longname 	= 'Fresh water transport (Antarctic Bottom Water)'
    
    fh.variables['time'].units 		    = 'Year'
    fh.variables['Transport'].units 	= 'Sv'
    fh.variables['F_OV'].units 		    = 'Sv'
    fh.variables['F_OV_ASW'].units 		= 'Sv'
    fh.variables['F_OV_AIW'].units 		= 'Sv'
    fh.variables['F_OV_NADW'].units 	= 'Sv'
    fh.variables['F_OV_ABW'].units 		= 'Sv'
    
    #Writing data to correct variable	
    fh.variables['time'][:]     	  	= time
    fh.variables['Transport'][:]    	= transport_all
    fh.variables['F_OV'][:] 			= transport_salt_all
    fh.variables['F_OV_ASW'][:] 		= transport_salt_ASW_all
    fh.variables['F_OV_AIW'][:] 		= transport_salt_AIW_all
    fh.variables['F_OV_NADW'][:] 		= transport_salt_NADW_all
    fh.variables['F_OV_ABW'][:] 		= transport_salt_ABW_all
    
    fh.close()
