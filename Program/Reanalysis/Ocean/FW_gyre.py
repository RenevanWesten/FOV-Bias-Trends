#Program determines the azonal (gyre) contribution

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

section_name	= 'FOV_section_34S'
#section_name	= 'FOV_section_60N'

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
layer_field_area  = ma.masked_all(shape(layer_field))

for depth_i in range(len(depth)):
	#Determine the surface area
	layer_field_area[depth_i]  	= layer_field[depth_i] * grid_x
   
#-----------------------------------------------------------------------------------------

#Define empty array's
transport_gyre_all	= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon, lat, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[file_i], lat_index, depth_min_index, depth_max_index, section_name)

	#Determine the zonal means
	v_vel_zonal = np.mean(v_vel, axis = 1)
	salt_zonal	= np.mean(salt, axis = 1)
	
	v_vel_prime	= ma.masked_all(np.shape(v_vel))
	salt_prime	= ma.masked_all(np.shape(salt))
	
	for depth_i in range(len(depth)):
		#Determine the differences with respect to the zonal means
		v_vel_prime[depth_i]	= v_vel[depth_i] - v_vel_zonal[depth_i]
		salt_prime[depth_i]	= salt[depth_i] - salt_zonal[depth_i]

	#Now determine the azonal component (gyre, in Sv)
	transport_gyre_all[file_i]	= (-1.0 / 35.0) * np.sum(v_vel_prime * salt_prime * layer_field_area) / 10**6.0

#-----------------------------------------------------------------------------------------

plot(time, transport_gyre_all, '-k')

show()

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/FW_gyre_'+section_name[4:]+'.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('F_gyre', float, ('time'), zlib=True)

fh.variables['F_gyre'].long_name 	= 'Freshwater transport by gyre'

fh.variables['time'].units 		= 'Year'
fh.variables['F_gyre'].units 		= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	= time
fh.variables['F_gyre'][:] 		= transport_gyre_all

fh.close()
