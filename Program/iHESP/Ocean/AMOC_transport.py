#Program determines the AMOC strength at 26N

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

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon 	= fh.variables['ULONG'][:]							                #Longitude
	lat 	= fh.variables['ULAT'][lat_index]						            #Latitude 
	depth   = fh.variables['z_t'][depth_min_index:depth_max_index] 			#Depth (m)
	layer	= fh.variables['dz'][depth_min_index:depth_max_index] 				#Layer thickness (m)
	grid_x	= fh.variables['DXU'][lat_index] 		                            #Zonal grid cell length (m)
	v_vel 	= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index] #Meridional velocity (m/s)

	fh.close()
    
	return lon, lat, depth, layer, grid_x, v_vel
    			
#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 1000

experiment  	= 'PI_control_HR'
#experiment  	= 'RCP_HR'
#experiment  	= 'PI_control_LR'
#experiment  	= 'RCP_LR'

#-----------------------------------------------------------------------------------------

files = glob.glob(directory+'Data/'+experiment+'/AMOC_section_26N/iHESP_data_year_*.nc')
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

lat 		= fh.variables['ULAT'][:]	#Latitude  
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lat_index		    = 1
depth_u		        = depth_u[lat_index]

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon, lat, depth, layer, grid_x, v_vel	= ReadinData(files[0], lat_index, depth_min_index, depth_max_index)
layer_field				= ma.masked_all((len(depth), len(lon)))
depth_top				= np.zeros(len(depth))

for depth_i in range(1, len(depth)):
	#Determine the depth top
	depth_top[depth_i]	= depth_top[depth_i-1] + layer[depth_i-1]

if experiment == 'PI_control_LR' or experiment == 'RCP_LR':
	layer_field	= ma.masked_all((len(depth), len(lon[0])))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field[depth_i]	= layer[depth_i]
	layer_field[depth_i]	= ma.masked_array(layer_field[depth_i], mask = v_vel[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff	= np.sum(layer_field, axis = 0) - depth_u

	if depth_i == len(depth) - 1:
		#Last layer, get the depth difference with respect to top and depth max boundary
		depth_diff	= layer_field[depth_i] - (depth_max -  depth_top[depth_i])

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff	= ma.masked_where(depth_diff < 0, depth_diff)
	depth_diff	= depth_diff.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field[depth_i]	= layer_field[depth_i] - depth_diff

#-----------------------------------------------------------------------------------------
transport_all		= ma.masked_all(len(time))

for file_i in range(len(files)):
	#Now determine for AMOC strength
	print(file_i)
	    
	lon, lat, depth, layer, grid_x, v_vel = ReadinData(files[file_i], lat_index, depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport	= v_vel * layer_field * grid_x

	#Determine the transport per depth layer (in Sv) and take sum to determine total transport
	transport_all[file_i]	= np.sum(transport) / 1000000.0
		        
#-----------------------------------------------------------------------------------------


fig, ax	= subplots()

graph_control	= ax.plot(time, transport_all, '-k', linewidth = 1.5)

ax.set_xlabel('Model year')
ax.set_ylabel('Volume transport (Sv)')
ax.set_ylim(0, 25)
ax.grid()

show()

#-----------------------------------------------------------------------------------------
#Do not overwrite the data, the data in the example directory is not the complete time series
#The time series contains all the availble model years
sys.exit()	

print('Data is written to file')
fh = netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_'+experiment+'.nc', 'w')

fh.createDimension('time', len(time))

fh.createVariable('time', float, ('time'), zlib=True)
fh.createVariable('Transport', float, ('time'), zlib=True)

fh.variables['Transport'].long_name 	= 'Volume transport'

fh.variables['time'].units 		    = 'Year'
fh.variables['Transport'].units 	= 'Sv'

#Writing data to correct variable	
fh.variables['time'][:]     	  	= time
fh.variables['Transport'][:]    	= transport_all

fh.close()
