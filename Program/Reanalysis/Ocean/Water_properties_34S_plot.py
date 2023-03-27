#Program determines the MOV index

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	   = '../../../Data/Reanalysis/'

def ReadinData(filename, lat_index, depth_min_index, depth_max_index, section_name):

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon 		= fh.variables['longitude'][:]						#Longitude
	lat 		= fh.variables['latitude'][lat_index]					#Latitude 
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

lat_FOV 	= -34

#-----------------------------------------------------------------------------------------

if lat_FOV == -34:
	section_name	= 'FOV_section_34S'

if lat_FOV == 26:
	section_name	= 'AMOC_section_26N'

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
	grid_x_depth          	= ma.masked_array(grid_x, mask = v_vel[depth_i].mask)
	grid_x_norm[depth_i]  	= grid_x_depth / np.sum(grid_x_depth)

#-----------------------------------------------------------------------------------------

#Define empty array's
vel_all			= ma.masked_all((len(time), len(depth)))
vel_salt_all    = ma.masked_all((len(time), len(depth)))
salt_all		= ma.masked_all((len(time), len(depth), len(lon)))

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

	#Save the meridional baroclinic transport
	vel_all[file_i]		= np.sum(vel_baroclinic * grid_x_norm, axis = 1) * 100.0
	vel_salt_all[file_i]	= (-1.0 / 35.0) * transport_clin * salt_zonal / 10**6.0
	salt_all[file_i]	= salt

vel_all			= np.mean(vel_all, axis = 0)
vel_salt_all		= np.mean(vel_salt_all, axis = 0)
vel_salt_all		= vel_salt_all / layer_field[:, 780] * 1000.0
salt_all		= np.mean(salt_all, axis = 0)

#-----------------------------------------------------------------------------------------
#Get the water properties

#North Atlantic Deep Water (NADW) has negative meridional velocities
depth_index_NADW = np.where((depth >= 700) & (vel_all <= 0))[0][0]

#Antarctic bottom water (ABW) is directly below the NADW, get the first index
depth_index_ABW	= np.where((depth >= 3000) & (vel_all >= 0))[0][0]

#The Antarctic Intermediate water is between the NADW and 500 m
depth_index_AIW	= np.where(depth >= 500)[0][0]


depth_top	= np.zeros(len(depth))

for depth_i in range(1, len(depth)):
	depth_top[depth_i]	= depth_top[depth_i - 1] + layer_field[depth_i - 1, 780]

depth_AIW	= depth_top[depth_index_AIW]
depth_NADW	= depth_top[depth_index_NADW]
depth_ABW	= depth_top[depth_index_ABW]

lon_AIW_index		= np.where(salt_all[depth_index_AIW].mask == False)[0]
lon_NADW_index		= np.where(salt_all[depth_index_NADW].mask == False)[0]
lon_ABW_index		= np.where(salt_all[depth_index_ABW].mask == False)[0]
lon_AIW_1, lon_AIW_2	= lon[lon_AIW_index[0]], lon[lon_AIW_index[-1]]
lon_NADW_1, lon_NADW_2	= lon[lon_NADW_index[0]], lon[lon_NADW_index[-1]]
lon_ABW_1, lon_ABW_2	= lon[lon_ABW_index[0]], lon[lon_ABW_index[-1]]

#-----------------------------------------------------------------------------------------

depth_crop			        = 1000
factor_depth_crop		    = 4
depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

if depth_AIW > depth_crop:
	depth_AIW	= ((depth_AIW - depth_crop) / factor_depth_crop) + depth_crop
if depth_NADW > depth_crop:
	depth_NADW	= ((depth_NADW - depth_crop) / factor_depth_crop) + depth_crop
if depth_ABW > depth_crop:
	depth_ABW	= ((depth_ABW - depth_crop) / factor_depth_crop) + depth_crop
#-----------------------------------------------------------------------------------------

if True:
    #For the CMIP6 plot
	fig, ax	= subplots()

	ax.fill_between([-60, 20], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)
	ax.plot([lon_AIW_1, lon_AIW_2], [depth_AIW, depth_AIW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_NADW_1, lon_NADW_2], [depth_NADW, depth_NADW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_ABW_1, lon_ABW_2], [depth_ABW, depth_ABW], linestyle = '--', linewidth = 2.0, color = 'k')

	CS	= contourf(lon, depth, salt_all, levels = np.arange(34, 36.01, 0.1), extend = 'both', cmap = 'BrBG_r')
	cbar	= colorbar(CS, ticks = np.arange(34, 36.01, 0.5))
	cbar.set_label('Salinity (g kg$^{-1}$)')

	ax.set_xlim(-60, 20)
	ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
	ax.set_ylabel('Depth (m)')	

	ax.set_xticks(np.arange(-60, 21, 10))
	ax.set_xticklabels(['60$^{\circ}$W', '50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

	labels =  ax.get_yticks()
	for label_i in range(len(labels)):
		if labels[label_i] > depth_crop:
			#Rescale the xlabels
			labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

	labels	= labels.astype(int)
	ax.set_yticklabels(labels)

	#Now get FOV and AMOC
	fh 	= netcdf.Dataset(directory+'Ocean/FOV_index_section_34S.nc')
	FOV	= np.mean(fh.variables['F_OV'][:27])		#FOV
	fh.close()


	fh 	= netcdf.Dataset(directory+'Ocean/AMOC_transport_depth_0-1000m.nc')
	AMOC	= np.mean(fh.variables['Transport'][:27])	#AMOC strength
	fh.close()

	ax.set_title('Reanalysis, $F_{\mathrm{ovS}}$ = '+str(round(FOV, 2))+' Sv, AMOC = '+str(round(AMOC, 1))+' Sv')

	#-----------------------------------------------------------------------------------------
	cNorm  		= colors.Normalize(vmin=34, vmax= 36) 		#Probablility
	scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='BrBG_r') 	#Using colormap
	color_fresh 	= scalarMap.to_rgba(34.5)
	color_salt 	= scalarMap.to_rgba(35.5)

	ax2 = fig.add_axes([0.125, 0.11, 0.05, 0.768])

	ax2.plot(vel_salt_all, depth, '-k', linewidth = 2.0)
	ax2.set_xlim(-0.1, 0.1)
	ax2.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)

	ax2.axvline(x = 0, linestyle = '--', color = 'k')
	ax2.fill_betweenx(depth, vel_salt_all, where = vel_salt_all >= 0.0, color = color_fresh, alpha = 0.50)	
	ax2.fill_betweenx(depth, vel_salt_all, where = vel_salt_all <= 0.0, color = color_salt, alpha = 0.50)

	labels =  ax.get_yticks()

	for label_i in labels:
		ax2.axhline(y = label_i, color = 'gray', linestyle = ':', alpha = 0.5)
		
	ax2.set_xticklabels([])
	ax2.set_yticklabels([])

	show()

#-----------------------------------------------------------------------------------------

cNorm  		= colors.Normalize(vmin=-1, vmax= 1) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='RdBu_r') 	#Using colormap
color_south 	= scalarMap.to_rgba(-0.5)
color_north 	= scalarMap.to_rgba(0.5)

fig, ax	= subplots()

ax.axhline(y = depth_AIW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.axhline(y = depth_NADW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.axhline(y = depth_ABW, linestyle = '--', linewidth = 2.0, color = 'k')
plot(vel_all, depth, '-k', linewidth = 2.0)

ax.set_xlim(-2, 2)
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

ax.set_xlabel('Meridional velocity (cm s$^{-1}$)')
ax.set_ylabel('Depth (m)')
ax.axvline(x = 0, linestyle = '--', color = 'k')

ax.text(1.9, 350, 'ASW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.9, 850, 'AIW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.9, 1350, 'NADW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.9, 1900, 'ABW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)

ax.set_title('a) Meridional velocity, Reanalysis (1994 - 2020)')

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()


ax.fill_between([-60, 20], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)
ax.plot([lon_AIW_1, lon_AIW_2], [depth_AIW, depth_AIW], linestyle = '--', linewidth = 2.0, color = 'k')
ax.plot([lon_NADW_1, lon_NADW_2], [depth_NADW, depth_NADW], linestyle = '--', linewidth = 2.0, color = 'k')
ax.plot([lon_ABW_1, lon_ABW_2], [depth_ABW, depth_ABW], linestyle = '--', linewidth = 2.0, color = 'k')

CS	= contourf(lon, depth, salt_all, levels = np.arange(34, 36.01, 0.1), extend = 'both', cmap = 'BrBG_r')
cbar	= colorbar(CS, ticks = np.arange(34, 36.01, 0.5))
cbar.set_label('Salinity (g kg$^{-1}$)')

ax.set_xlim(-60, 20)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.set_ylabel('Depth (m)')	

ax.set_xticks(np.arange(-60, 21, 10))
ax.set_xticklabels(['60$^{\circ}$W', '50$^{\circ}$W', '40$^{\circ}$W', '30$^{\circ}$W', '20$^{\circ}$W', '10$^{\circ}$W','0$^{\circ}$', '10$^{\circ}$E', '20$^{\circ}$E'])

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)


ax.text(-18, 350, 'ASW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
ax.text(-18, 850, 'AIW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
ax.text(-18, 1350, 'NADW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)
ax.text(-18, 1900, 'ABW', verticalalignment='center', horizontalalignment='center', color = 'k', fontsize=16)


ax.set_title('d) Salinity, Reanalysis (1994 - 2020)')

#-----------------------------------------------------------------------------------------

cNorm  		= colors.Normalize(vmin=34, vmax= 36) 		#Probablility
scalarMap 	= cm.ScalarMappable(norm=cNorm, cmap='BrBG_r') 	#Using colormap
color_fresh 	= scalarMap.to_rgba(34.5)
color_salt 	= scalarMap.to_rgba(35.5)

fig, ax	= subplots()

ax.axhline(y = depth_AIW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.axhline(y = depth_NADW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.axhline(y = depth_ABW, linestyle = '--', linewidth = 2.0, color = 'k')
ax.plot(vel_salt_all, depth, '-k', linewidth = 2.0)

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(((5500 - depth_crop) / factor_depth_crop) + depth_crop, 0)
ax.grid()

labels =  ax.get_yticks()
for label_i in range(len(labels)):
	if labels[label_i] > depth_crop:
		#Rescale the xlabels
		labels[label_i]	= ((labels[label_i] - depth_crop) * factor_depth_crop) + depth_crop

labels	= labels.astype(int)
ax.set_yticklabels(labels)

ax.set_xlabel(r'Freshwater transport (mSv m$^{-1}$)')
ax.set_ylabel('Depth (m)')
ax.axvline(x = 0, linestyle = '--', color = 'k')

ax.fill_betweenx(depth, vel_salt_all, where = vel_salt_all >= 0.0, color = color_fresh, alpha = 0.50)	
ax.fill_betweenx(depth, vel_salt_all, where = vel_salt_all <= 0.0, color = color_salt, alpha = 0.50)

ax.text(1.45, 350, 'ASW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.45, 850, 'AIW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.45, 1350, 'NADW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)
ax.text(1.45, 1900, 'ABW', verticalalignment='center', horizontalalignment='right', color = 'k', fontsize=16)	

ax.set_title('g) Freshwater transport, Reanalysis (1994 - 2020)')

show()
#-----------------------------------------------------------------------------------------

