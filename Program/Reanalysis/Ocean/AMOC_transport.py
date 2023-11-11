#Program determines the AMOC strength at 26N

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
depth_max	= 1000

section_name	= 'AMOC_section_26N'

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

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lat 		= fh.variables['latitude'][:]	#Latitude  
depth   	= fh.variables['depth'][:]	#Depth (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lat_index		    = (fabs(26 - lat)).argmin()

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, depth, layer_field, grid_x, v_vel, salt = ReadinData(files[0], lat_index, depth_min_index, depth_max_index, section_name)

for lon_i in range(len(lon)):
	#Get all the layers which have a maximum depth below given range
	if np.sum(layer_field[:, lon_i]) > depth_max:
		#Adjust the last layer
		layer_field[-1, lon_i]	-= (np.sum(layer_field[:, lon_i]) - depth_max)

#-----------------------------------------------------------------------------------------

#Define empty array's
transport_all		= ma.masked_all(len(time))

for time_i in range(len(time)):
	#Now determine for each month
	print(time_i)
	    
	lon, lat, depth, layer_field_old, grid_x, v_vel, salt = ReadinData(files[time_i], lat_index, depth_min_index, depth_max_index, section_name)

	#Determine the meridional transport
	transport		= v_vel * layer_field * grid_x

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[time_i]	= np.sum(transport) / 1000000.0

#-----------------------------------------------------------------------------------------

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].longname 	= 'Volume transport'

fh.variables['time'].units 		= 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all

fh.close()
