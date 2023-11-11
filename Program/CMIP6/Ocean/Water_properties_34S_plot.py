#Program determines the water properties at 34S

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory	= '../../../Data/CMIP6/Ocean/'

def ReadinData(filename, lat_index, depth_min_index, depth_max_index):
	#The CESM grid is structured as
	#S - S -
	#- U* - U = 34S
	#S* - S - 
	#Where the stars have the same index

	fh = netcdf.Dataset(filename, 'r')

	#First get the u-grid
	lon_u 		= fh.variables['ULONG'][:]							        #Longitude
	lat_u 		= fh.variables['ULAT'][lat_index]							#Latitude 
	depth   	= fh.variables['z_t'][depth_min_index:depth_max_index] 				   	#Depth (m)
	layer		= fh.variables['dz'][depth_min_index:depth_max_index] 				   	#Layer thickness (m)
	grid_x_u	= fh.variables['DXU'][lat_index] 		                               		#Zonal grid cell length (m)
	v_vel 		= fh.variables['VVEL'][depth_min_index:depth_max_index, lat_index] 			#Meridional velocity (m/s)

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

#-----------------------------------------------------------------------------------------

#Get the model names and path
models = glob.glob(directory+'FOV_section_34S/*.nc')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][len(directory)+32:-3]

#-----------------------------------------------------------------------------------------
for model_i in range(len(models)):
	#For each model get the all the files
	print('-----------------------------------------------------------------------------------------')
	print(str(model_i+1)+') '+models[model_i])
	print('-----------------------------------------------------------------------------------------')

	fh = netcdf.Dataset(directory+'FOV_section_34S/FOV_section_34S_'+models[model_i]+'.nc', 'r')

	if models[model_i] == 'CAS-ESM2-0' or models[model_i] == 'CIESM' or models[model_i] == 'GFDL-CM4':
		#These models have the same lon/lat-grid for the velocity and salinity
		time		= fh.variables['time'][:]
		lon		    = fh.variables['lon'][:] 	#longitude
		depth   	= fh.variables['depth'][:] 	#Depth (m)
		layer		= fh.variables['LAYER'][:] 	#Layer thickness (m)
		grid_x		= fh.variables['DX'][:]		#Zonal grid cell length (m)
		v_vel 		= fh.variables['VVEL'][:] 	#Meridional velocity (m/s)
		salt		= fh.variables['SALT'][:] 	#Salinity (g / kg)

		fh.close()

		#-----------------------------------------------------------------------------------------
		#Normalise layer field per layer
		grid_x_norm	= ma.masked_all((len(depth), len(lon)))

		for depth_i in range(len(depth)):
			#Normalise the length per depth level
			grid_x_depth          = ma.masked_array(grid_x, mask = salt[0, depth_i].mask)
			grid_x_norm[depth_i]  = grid_x_depth / np.sum(grid_x_depth)

		#-----------------------------------------------------------------------------------------

		#Define empty array's
		vel_all			= ma.masked_all((len(time), len(depth)))
		vel_salt_all		= ma.masked_all((len(time), len(depth)))
		salt_all		= ma.masked_all((len(time), len(depth), len(lon)))

		for year_i in range(len(time)):
			#Loop over each year

			#Determine the meridional transport
			transport	= v_vel[year_i] * layer * grid_x

			#Determine the section averaged velocity (barotropic)
			vel_barotropic	= np.sum(transport) / np.sum(layer * grid_x)

			#Determine the overturning velocity (baroclinic)
			vel_baroclinic	 = v_vel[year_i] - vel_barotropic

			#Determine the zonal means
			salt_zonal      = np.sum(salt[year_i] * grid_x_norm, axis = 1)  - 35.0
			transport_clin	= np.sum(vel_baroclinic * layer * grid_x, axis = 1)

			#Save the meridional baroclinic transport
			vel_all[year_i]			= np.sum(vel_baroclinic * grid_x_norm, axis = 1) * 100.0
			vel_salt_all[year_i]		= (-1.0 / 35.0) * transport_clin * salt_zonal / 10**6.0
			salt_all[year_i]		= salt[year_i]

		vel_all		= np.mean(vel_all, axis = 0)
		vel_salt_all	= np.mean(vel_salt_all, axis = 0) / np.max(layer, axis = 1) * 1000.0
		salt_all	= np.mean(salt_all, axis = 0)

		depth_top	= np.zeros(len(depth))

		for depth_i in range(1, len(depth)):
			depth_top[depth_i]	= depth_top[depth_i - 1] + np.max(layer[depth_i - 1])

	else:
		#First get the related variables for the velocity grid
		time		= fh.variables['time'][:]
		lon_vel		= fh.variables['lon_VEL'][:] 	#longitude
		depth   	= fh.variables['depth'][:] 	#Depth (m)
		layer_vel	= fh.variables['LAYER_VEL'][:] 	#Layer thickness (m)
		grid_x_vel	= fh.variables['DX_VEL'][:]	#Zonal grid cell length (m)
		v_vel 		= fh.variables['VVEL'][:] 	#Meridional velocity (m/s)

		#Next get the related variables for the salinity grid
		lon_salt	= fh.variables['lon_SALT'][:] 	#longitude
		layer_salt	= fh.variables['LAYER_SALT'][:] #Layer thickness (m)
		grid_x_salt	= fh.variables['DX_SALT'][:]	#Zonal grid cell length (m)
		salt		= fh.variables['SALT'][:] 	#Salinity (g / kg)

		fh.close()

		#-----------------------------------------------------------------------------------------
		#Normalise layer field per layer
		grid_x_vel_norm		= ma.masked_all((len(depth), len(lon_vel)))
		grid_x_salt_norm	= ma.masked_all((len(depth), 2, len(lon_salt)))
		lat_weight		= ma.masked_all((len(depth), 2))

		for depth_i in range(len(depth)):
			#Normalise the length per depth level
			grid_x_vel_depth          	= ma.masked_array(grid_x_vel, mask = v_vel[0, depth_i].mask)
			grid_x_salt_depth          	= ma.masked_array(grid_x_salt, mask = salt[0, depth_i].mask)

			#Normalise the vel depth
			grid_x_vel_norm[depth_i]  	= grid_x_vel_depth / np.sum(grid_x_vel_depth)

			#Now get the lat weights
			lat_weight[depth_i]		= np.sum(grid_x_salt_depth, axis = 1) / np.sum(grid_x_salt_depth)

			for lat_i in range(2):
				grid_x_salt_norm[depth_i, lat_i]  = grid_x_salt_depth[lat_i] / np.sum(grid_x_salt_depth[lat_i])


		#-----------------------------------------------------------------------------------------

		#Define empty array's
		vel_all		= ma.masked_all((len(time), len(depth)))
		vel_salt_all	= ma.masked_all((len(time), len(depth)))
		salt_all	= ma.masked_all((len(time), len(depth), len(lon_salt)))

		for year_i in range(len(time)):
			#Loop over each year

			#Determine the meridional transport
			transport	= v_vel[year_i] * layer_vel * grid_x_vel

			#Determine the section averaged velocity (barotropic)
			vel_barotropic	= np.sum(transport) / np.sum(layer_vel * grid_x_vel)

			#Determine the overturning velocity (baroclinic)
			vel_baroclinic	 = v_vel[year_i] - vel_barotropic

			#Determine the zonal means
			salt_zonal      = np.sum(salt[year_i] * grid_x_salt_norm, axis = 2)  - 35.0
			transport_clin	= np.sum(vel_baroclinic * layer_vel * grid_x_vel, axis = 1)

			#Take the mean over the two latitudes for the salt
			salt_zonal      = np.sum(salt_zonal * lat_weight, axis = 1)

			#Save the meridional baroclinic transport
			vel_all[year_i]			= np.sum(vel_baroclinic * grid_x_vel_norm, axis = 1) * 100.0
			vel_salt_all[year_i]		= (-1.0 / 35.0) * transport_clin * salt_zonal / 10**6.0
			salt_all[year_i]		= salt[year_i, :, 0]

		vel_all		= np.mean(vel_all, axis = 0)
		vel_salt_all	= np.mean(vel_salt_all, axis = 0) / np.max(layer_vel, axis = 1) * 1000.0
		salt_all	= np.mean(salt_all, axis = 0)
		lon		= np.copy(lon_salt)

		depth_top	= np.zeros(len(depth))

		for depth_i in range(1, len(depth)):
			depth_top[depth_i]	= depth_top[depth_i - 1] + np.max(layer_vel[depth_i - 1])

	#-----------------------------------------------------------------------------------------
	#Get the water properties

	#North Atlantic Deep Water (NADW) has negative meridional velocities
	depth_index_NADW = np.where((depth >= 700) & (vel_all <= 0))[0][0]

	#Antarctic bottom water (AABW) is directly below the NADW, get the first index
	depth_index_AABW	= np.where((depth >= 3000) & (vel_all >= 0))[0][0]

	#The Antarctic Intermediate water is between the NADW and 500 m
	depth_index_AAIW	= np.where(depth >= 500)[0][0]

	depth_AAIW	= depth_top[depth_index_AAIW]
	depth_NADW	= depth_top[depth_index_NADW]
	depth_AABW	= depth_top[depth_index_AABW]

	lon_AAIW_index		= np.where(salt_all[depth_index_AAIW].mask == False)[0]
	lon_NADW_index		= np.where(salt_all[depth_index_NADW].mask == False)[0]
	lon_AABW_index		= np.where(salt_all[depth_index_AABW].mask == False)[0]
	lon_AAIW_1, lon_AAIW_2	= lon[lon_AAIW_index[0]], lon[lon_AAIW_index[-1]]
	lon_NADW_1, lon_NADW_2	= lon[lon_NADW_index[0]], lon[lon_NADW_index[-1]]
	lon_AABW_1, lon_AABW_2	= lon[lon_AABW_index[0]], lon[lon_AABW_index[-1]]

	#-----------------------------------------------------------------------------------------

	depth_crop			= 1000
	factor_depth_crop		= 4
	depth[depth > depth_crop] 	= ((depth[depth > depth_crop] - depth_crop) / factor_depth_crop) + depth_crop

	if depth_AAIW > depth_crop:
		depth_AAIW	= ((depth_AAIW - depth_crop) / factor_depth_crop) + depth_crop
	if depth_NADW > depth_crop:
		depth_NADW	= ((depth_NADW - depth_crop) / factor_depth_crop) + depth_crop
	if depth_AABW > depth_crop:
		depth_AABW	= ((depth_AABW - depth_crop) / factor_depth_crop) + depth_crop

	#-----------------------------------------------------------------------------------------

	fig, ax	= subplots()


	ax.fill_between([-60, 20], y1 = np.zeros(2) + depth[0], y2 = np.zeros(2) + 2*depth[-1], color = 'gray', alpha = 0.50)
	ax.plot([lon_AAIW_1, lon_AAIW_2], [depth_AAIW, depth_AAIW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_NADW_1, lon_NADW_2], [depth_NADW, depth_NADW], linestyle = '--', linewidth = 2.0, color = 'k')
	ax.plot([lon_AABW_1, lon_AABW_2], [depth_AABW, depth_AABW], linestyle = '--', linewidth = 2.0, color = 'k')

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
	fh 	= netcdf.Dataset(directory+'FOV_index_'+models[model_i]+'.nc')
	FOV	= np.mean(fh.variables['F_OV'][:27])		#FOV
	fh.close()


	fh 	= netcdf.Dataset(directory+'AMOC_transport_depth_0-1000m_'+models[model_i]+'.nc')
	AMOC	= np.mean(fh.variables['Transport'][:27])	#AMOC strength
	fh.close()

	ax.set_title(models[model_i]+', $F_{\mathrm{ovS}}$ = '+str(round(FOV, 2))+' Sv, AMOC = '+str(round(AMOC, 1))+' Sv')

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
