#Program plots the Indonesian Throughflow (ITF)

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors
from cartopy import crs as ccrs, feature as cfeature
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as mticker
from scipy import stats

#Making pathway to folder with all data
directory	  = '../../../Data/iHESP/'

def ReadinDataFOV(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	#Freshwater (Sv)

	fh.close()

	return time, FOV

def ReadinData(filename, lon_index, depth_min_index, depth_max_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon_u 		= fh.variables['ULONG'][lon_index]							                #Longitude
	lat_u 		= fh.variables['ULAT'][:]							                        #Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 				        #Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 				        #Layer thickness (m)
	grid_y_u	= fh.variables['DYU'][:, lon_index] 		                               	#Meridional grid cell length (m)
	u_vel 		= fh.variables['UVEL'][depth_min_index:depth_max_index, :, lon_index] 		#Meridional velocity (m/s)

	#Get the t-grid
	lon_t 		= fh.variables['TLONG'][lon_index:lon_index+2]						        #Longitude
	lat_t 		= fh.variables['TLAT'][:]								                    #Latitude 
	grid_y_t	= fh.variables['DYT'][:, lon_index:lon_index+2] 		                   	#Meridional grid cell length (m)
	salt		= fh.variables['SALT'][depth_min_index:depth_max_index, :, lon_index:lon_index+2] 	#Salinity (g / kg)

	fh.close()
    
	return lon_u, lat_u, lon_t, lat_t, depth, layer, grid_y_u, grid_y_t, u_vel, salt
    

def SignificantTrend(time, data):
	"""Finds whether trend is significant
	Returns the trend and if it significant (= 1)"""

	#Set time similar to Santer et al. (2000), time array from 1 till N
	#Statistical significance of trends and trend differences in layer-average atmospheric salterature time series
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


#Read in the time series of the PI control and forcing (historical + RCP8.5)
time_control, FOV_control	= ReadinDataFOV(directory+'Ocean/ITF_PI_control_HR.nc')
time_rcp, FOV_rcp		    = ReadinDataFOV(directory+'Ocean/ITF_RCP_HR.nc')
time_control			    += 1599

print('Trend model year:', int(time_rcp[150]),'-', int(time_rcp[250]))
print()
trend, error, p_value = SignificantTrend(time_rcp[150:251], FOV_rcp[150:251])
print('ITF:', trend*100, '('+str(p_value)+')')
print()

#-----------------------------------------------------------------------------------------
depth_min 	= 0
depth_max	= 100

fh = netcdf.Dataset(directory+'Ocean/SALT_SO_'+str(depth_min)+'_'+str(depth_max)+'m_trend_year_21-100_PI_control_HR.nc', 'r')

lon_ocn		= fh.variables['lon'][:] 	
lat_ocn		= fh.variables['lat'][:] 		
salt_trend	= fh.variables['SALT_trend'][:] 	

fh.close()
#-----------------------------------------------------------------------------------------

depth_min 	= 0
depth_max	= 6000

files = glob.glob(directory+'Data/PI_control_HR/ITF/iHESP_data_year_0050.nc')
files.sort()

#-----------------------------------------------------------------------------------------

#Get all the relevant indices to determine the mass transport
fh = netcdf.Dataset(files[0], 'r')

lon_u 		= fh.variables['ULONG'][:]	#Longitude
depth   	= fh.variables['z_t'][:]	#Depth (m)
depth_u 	= fh.variables['HU'][:] 	#Depth at u-grid (m)
	
fh.close()

#Get the dimensions of depth and latitude
depth_min_index 	= (fabs(depth_min - depth)).argmin()
depth_max_index 	= (fabs(depth_max - depth)).argmin() + 1
lon_index		    = 1
depth_u		        = depth_u[:, lon_index]

#-----------------------------------------------------------------------------------------
#Determine the section length per depth layer
lon_u, lat_u, lon_t, lat_t, depth, layer, grid_y_u, grid_y_t, u_vel, salt 	= ReadinData(files[0], lon_index, depth_min_index, depth_max_index)
layer_field_u	= ma.masked_all((len(depth), len(lat_u)))

for depth_i in range(len(depth)):
	#Determine the total length of the section, based on non-masked elements
	layer_field_u[depth_i]	= layer[depth_i]
	layer_field_u[depth_i]	= ma.masked_array(layer_field_u[depth_i], mask = u_vel[depth_i].mask)

	#Determine where the layer needs to be adjusted, partial depth cells
	depth_diff_u	= np.sum(layer_field_u, axis = 0) - depth_u

	#If the depth difference is negative (i.e. bottom is not reached), set to zero
	depth_diff_u	= ma.masked_where(depth_diff_u < 0, depth_diff_u)
	depth_diff_u	= depth_diff_u.filled(fill_value = 0.0)

	#Subtract the difference of the current layer with the difference
	layer_field_u[depth_i]	= layer_field_u[depth_i] - depth_diff_u

#Normalise layer field per layer
layer_field_u_norm  = ma.masked_all(shape(layer_field_u))
grid_y_u_norm	    = ma.masked_all((len(depth), len(lat_u)))
grid_y_t_norm	    = ma.masked_all((len(depth), len(lat_t), 2))
lon_weight	        = ma.masked_all((len(depth), 2))

for depth_i in range(len(depth)):
	#Normalise each layer
	layer_field_u_norm[depth_i]		= layer_field_u[depth_i] / np.sum(layer_field_u[depth_i])
    
	#Normalise the length
	grid_y_u_depth          		= ma.masked_array(grid_y_u, mask = u_vel[depth_i].mask)
	grid_y_u_norm[depth_i]  		= grid_y_u_depth / np.sum(grid_y_u_depth)
	grid_y_t_depth          		= ma.masked_array(grid_y_t, mask = salt[depth_i].mask)

	#Now get the lon weights
	lon_weight[depth_i]			= np.sum(grid_y_t_depth, axis = 0) / np.sum(grid_y_t_depth)

	for lon_i in range(2):
		grid_y_t_norm[depth_i, :, lon_i]  = grid_y_t_depth[:, lon_i] / np.sum(grid_y_t_depth[:, lon_i])

#-----------------------------------------------------------------------------------------

#Define empty array's
vel_all		= ma.masked_all((len(files), len(depth), 1))
salt_all	= ma.masked_all((len(files), len(depth), len(lat_t)))

for file_i in range(len(files)):
	#Now determine for each month
	print(file_i)
	    
	lon_u, lat_u, lon_t, lat_t, depth, layer, grid_y_u, grid_y_t, u_vel, salt = ReadinData(files[file_i], lon_index, depth_min_index, depth_max_index)

	#Determine the meridional transport
	transport	= u_vel * layer_field_u * grid_y_u

	#Determine the section averaged velocity (barotropic)
	vel_barotropic	= np.sum(transport) / np.sum(layer_field_u * grid_y_u)

	#Determine the overturning velocity (baroclinic)
	vel_baroclinic	 = u_vel - vel_barotropic

	#Save the meridional baroclinic transport
	vel_all[file_i, :, 0]		= np.sum(vel_baroclinic * grid_y_u_norm, axis = 1) * 100.0
	salt_all[file_i]		= salt[:, :, 0]

vel_all			= np.mean(vel_all[:, :, 0], axis = 0)
salt_all		= np.mean(salt_all, axis = 0)

#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

graph_control	= ax.plot(time_control, FOV_control, '-k', linewidth = 1.5, label = 'PI control')
graph_rcp	= ax.plot(time_rcp, FOV_rcp, '-r', linewidth = 1.5, label = 'Hist/RCP8.5')

ax.set_xlabel('Model year')
ax.set_ylabel('Freshwater transport (Sv)')
ax.set_xlim(1600, 2100)
ax.set_ylim(-0.3, 0.3)
ax.grid()
ax.set_xticks([1600, 1700, 1800, 1900, 2000, 2100])
ax.set_xticklabels(['1', '100', '200', '300/1900', '400/2000',  '500/2100'])

graphs	      = graph_control + graph_rcp

legend_labels = [l.get_label() for l in graphs]
ax.legend(graphs, legend_labels, loc='lower left', ncol=1, framealpha = 1.0, numpoints = 1)

ax.set_title('a) $F_{\mathrm{ov}}$ (into the Indian Ocean) at 115$^{\circ}$E, HR-CESM')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.35, 0.18, 0.50, 0.30], projection = ccrs.PlateCarree())

CS      = ax2.contourf(lon_ocn, lat_ocn, salt_trend, levels = np.arange(-1, 1.01, 0.05), extend = 'both', cmap = 'BrBG_r', transform=ccrs.PlateCarree())

cbar    = colorbar(CS, ticks = np.arange(-1, 1.01, 1), fraction=0.021, pad=0.04)
cbar.set_label('Salinity trend (g kg$^{-1}$ per century)', fontsize = 8)

gl = ax2.gridlines(draw_labels=True)
gl.top_labels = False
gl.right_labels = False
gl.xlocator = mticker.FixedLocator([100, 110, 120, 130, 140])
gl.ylocator = mticker.FixedLocator([-20, -10, 0])
ax2.set_extent([99.9999, 140.0001, -25, 0.0001], ccrs.PlateCarree())
ax2.coastlines('50m')
ax2.add_feature(cfeature.LAND, zorder=0)
ax2.set_title('Salinity (0 - 100 m) trend, PI control (21 - 100)', fontsize = 10)

y_1	= np.arange(-23, -7.99, 0.1)
x_1	= np.zeros(len(y_1)) + 115
x_2	= np.arange(113.5, 116.51, 0.1)
y_2	= np.zeros(len(x_2)) + y_1[0]
x_3	= np.arange(113.5, 116.51, 0.1)
y_3	= np.zeros(len(x_3)) + y_1[-1]

ax2.plot(x_1, y_1, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_2, y_2, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)
ax2.plot(x_3, y_3, '-k', linewidth = 2.0, transform=ccrs.PlateCarree(), zorder = 10)

#-----------------------------------------------------------------------------------------

cNorm  		= colors.Normalize(vmin=-1, vmax= 1) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='RdBu_r') 	#Using colormap
color_south 	= scalarMap.to_rgba(-0.5)
color_north 	= scalarMap.to_rgba(0.5)

fig, ax	= subplots()

ax.plot(vel_all, depth, '-k', linewidth = 2.0)

ax.set_xlim(-6, 6)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.grid()

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.fill_betweenx(depth, vel_all, where = vel_all >= 0.0, color = color_north, alpha = 0.50)	
ax.fill_betweenx(depth, vel_all, where = vel_all <= 0.0, color = color_south, alpha = 0.50)	

ax.set_xlabel('Zonal velocity (cm s$^{-1}$)')
ax.set_ylabel('Depth (m)')
ax.axvline(x = 0, linestyle = '--', color = 'k')

ax.set_title('c) Zonal velocity, HR-CESM, PI control (model year 50)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-23, -8], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)

CS	= contourf(lat_t, depth, salt_all, levels = np.arange(34, 36.01, 0.1), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(34, 36.01, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-23, -8)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-20, -9, 5))
ax.set_xticklabels(['20$^{\circ}$S', '15$^{\circ}$S', '10$^{\circ}$S'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_title('e) Salinity, HR-CESM, PI control (model year 50)')

show()
#-----------------------------------------------------------------------------------------

