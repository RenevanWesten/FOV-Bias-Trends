#Program plots the freshwater transport and AMOC strength in CMIP6

from pylab import *
import numpy
import datetime
import time
import glob, os
import math
import netCDF4 as netcdf
import matplotlib.colors as colors

#Making pathway to folder with all data
directory_CMIP6	    = '../../../Data/CMIP6/Ocean/'
directory_iHESP     = '../../../Data/iHESP/Ocean/'
directory_reanalysis= '../../../Data/Reanalysis/Ocean/'

def ReadinData(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	FOV		    = fh.variables['F_OV'][:]	    #Freshwater (Sv)
	FOV_ASW		= fh.variables['F_OV_ASW'][:]	#Freshwater (Sv)
	FOV_AIW		= fh.variables['F_OV_AIW'][:]	#Freshwater (Sv)
	FOV_NADW	= fh.variables['F_OV_NADW'][:]	#Freshwater (Sv)
	FOV_ABW		= fh.variables['F_OV_ABW'][:]	#Freshwater (Sv)

	fh.close()

	return time, FOV, FOV_ASW, FOV_AIW, FOV_NADW, FOV_ABW

def ReadinDataAMOC(filename):

	fh = netcdf.Dataset(filename, 'r')

	time		= fh.variables['time'][:]		
	AMOC		= fh.variables['Transport'][:]	#AMOC strength

	fh.close()

	return time, AMOC

#-----------------------------------------------------------------------------------------
#--------------------------------MAIN SCRIPT STARTS HERE----------------------------------
#-----------------------------------------------------------------------------------------	

depth_min 	= 0
depth_max	= 1000

#-----------------------------------------------------------------------------------------

#Get the model names and path
models = glob.glob(directory_CMIP6+'FOV_index_*')
models.sort()

for model_i in range(len(models)):
	#Only retain the model names
	models[model_i]	= models[model_i][len(directory_CMIP6)+10:-3]

transport_salt_all	= ma.masked_all((len(models), 2))
transport_salt_ASW_all	= ma.masked_all((len(models), 2))
transport_salt_AIW_all	= ma.masked_all((len(models), 2))
transport_salt_NADW_all	= ma.masked_all((len(models), 2))
transport_salt_ABW_all	= ma.masked_all((len(models), 2))
AMOC_all		= ma.masked_all((len(models), 2))

transport_salt_all_time	= ma.masked_all((107, len(models)))
transport_salt_ASW_time	= ma.masked_all((107, len(models)))
transport_salt_AIW_time	= ma.masked_all((107, len(models)))
transport_salt_NADW_time= ma.masked_all((107, len(models)))
AMOC_all_time		= ma.masked_all((107, len(models)))

for model_i in range(len(models)):
	#Now read-in all the data for each CMIP6 model
	filename_FOV 					= directory_CMIP6+'FOV_index_'+models[model_i]+'.nc'
	filename_AMOC 					= directory_CMIP6+'AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_'+models[model_i]+'.nc'
	time, FOV, FOV_ASW, FOV_AIW, FOV_NADW, FOV_ABW	= ReadinData(filename_FOV)
	time, AMOC					= ReadinDataAMOC(filename_AMOC)

	#Save the time series
	transport_salt_all_time[:len(time), model_i]	= FOV
	transport_salt_ASW_time[:len(time), model_i]	= FOV_ASW
	transport_salt_AIW_time[:len(time), model_i]	= FOV_AIW
	transport_salt_NADW_time[:len(time), model_i]	= FOV_NADW
	AMOC_all_time[:len(time), model_i]		= AMOC

	#Regions over which we determine the time mean and trend
	time_present_start	= 1994
	time_present_end	= 2020
	time_future_start	= 2000
	time_future_end		= 2100

	time_present_start	= (np.abs(time - time_present_start)).argmin()
	time_present_end	= (np.abs(time - time_present_end)).argmin()+1
	time_future_start	= (np.abs(time - time_future_start)).argmin()
	time_future_end		= (np.abs(time - time_future_end)).argmin()+1

	#Determine the present-day time mean
	transport_salt_all[model_i, 0]		= np.mean(FOV[time_present_start:time_present_end])
	transport_salt_ASW_all[model_i, 0]	= np.mean(FOV_ASW[time_present_start:time_present_end])
	transport_salt_AIW_all[model_i, 0]	= np.mean(FOV_AIW[time_present_start:time_present_end])
	transport_salt_NADW_all[model_i, 0]	= np.mean(FOV_NADW[time_present_start:time_present_end])
	transport_salt_ABW_all[model_i, 0]	= np.mean(FOV_ABW[time_present_start:time_present_end])
	AMOC_all[model_i, 0]			= np.mean(AMOC[time_present_start:time_present_end])

	#Determine the trends
	FOV_trend, b		= np.polyfit(time[time_future_start:time_future_end], FOV[time_future_start:time_future_end], 1)
	FOV_ASW_trend, b	= np.polyfit(time[time_future_start:time_future_end], FOV_ASW[time_future_start:time_future_end], 1)
	FOV_AIW_trend, b	= np.polyfit(time[time_future_start:time_future_end], FOV_AIW[time_future_start:time_future_end], 1)
	FOV_NADW_trend, b	= np.polyfit(time[time_future_start:time_future_end], FOV_NADW[time_future_start:time_future_end], 1)
	FOV_ABW_trend, b	= np.polyfit(time[time_future_start:time_future_end], FOV_ABW[time_future_start:time_future_end], 1)
	AMOC_trend, b		= np.polyfit(time[time_future_start:time_future_end], AMOC[time_future_start:time_future_end], 1)

	#Save the trend (per century)
	transport_salt_all[model_i, 1]		= FOV_trend * 100.0
	transport_salt_ASW_all[model_i, 1]	= FOV_ASW_trend * 100.0
	transport_salt_AIW_all[model_i, 1]	= FOV_AIW_trend * 100.0
	transport_salt_NADW_all[model_i, 1]	= FOV_NADW_trend * 100.0
	transport_salt_ABW_all[model_i, 1]	= FOV_ABW_trend * 100.0
	AMOC_all[model_i, 1]			= AMOC_trend * 100.0

#Check the models with a correct FOV/AMOC value
models_1	= np.where((transport_salt_all[:, 0] >= -0.28) & (transport_salt_all[:, 0] <= -0.05))[0]
models_2	= np.where((AMOC_all[:, 0] >= 16) & (AMOC_all[:, 0] <= 19))[0]
models_3	= np.where((transport_salt_all[:, 0] > -0.05) & (np.abs(AMOC_all[:, 0] - 17.5)  > 1.5))[0]

print('Realistic FOV ('+str(len(models_1))+'):', np.mean(transport_salt_all[models_1, 0]), np.mean(AMOC_all[models_1, 0]))
print('Realistic AMOC ('+str(len(models_2))+'):', np.mean(transport_salt_all[models_2, 0]), np.mean(AMOC_all[models_2, 0]))
print('No realistic FOV/AMOC ('+str(len(models_3))+'):', np.mean(transport_salt_all[models_3, 0]), np.mean(AMOC_all[models_3, 0]))
print()

#-----------------------------------------------------------------------------------------
#Read in the time series of the HR-CESM forcing (historical + RCP8.5)
time_HR, FOV_HR, FOV_ASW_HR, FOV_AIW_HR, FOV_NADW_HR, FOV_ABW_HR = ReadinData(directory_iHESP+'FOV_index_section_34S_RCP_HR.nc')
time_HR, AMOC_HR						 = ReadinDataAMOC(directory_iHESP+'AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP_HR.nc')
time_LR, FOV_LR, FOV_ASW_LR, FOV_AIW_LR, FOV_NADW_LR, FOV_ABW_LR = ReadinData(directory_iHESP+'FOV_index_section_34S_RCP_LR.nc')
time_LR, AMOC_LR						 = ReadinDataAMOC(directory_iHESP+'AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m_RCP_LR.nc')

#Regions over which we determine the time mean and trend
time_present_start	= 1994
time_present_end	= 2020
time_future_start	= 2000
time_future_end		= 2100

time_present_start	= (np.abs(time_HR - time_present_start)).argmin()
time_present_end	= (np.abs(time_HR - time_present_end)).argmin()+1
time_future_start	= (np.abs(time_HR - time_future_start)).argmin()
time_future_end		= (np.abs(time_HR - time_future_end)).argmin()+1

#Determine the present-day time mean
FOV_mean_HR		= np.mean(FOV_HR[time_present_start:time_present_end])
FOV_ASW_mean_HR		= np.mean(FOV_ASW_HR[time_present_start:time_present_end])
FOV_AIW_mean_HR		= np.mean(FOV_AIW_HR[time_present_start:time_present_end])
FOV_NADW_mean_HR	= np.mean(FOV_NADW_HR[time_present_start:time_present_end])
FOV_ABW_mean_HR		= np.mean(FOV_ABW_HR[time_present_start:time_present_end])
AMOC_mean_HR		= np.mean(AMOC_HR[time_present_start:time_present_end])

FOV_mean_LR		= np.mean(FOV_LR[time_present_start:time_present_end])
FOV_ASW_mean_LR		= np.mean(FOV_ASW_LR[time_present_start:time_present_end])
FOV_AIW_mean_LR		= np.mean(FOV_AIW_LR[time_present_start:time_present_end])
FOV_NADW_mean_LR	= np.mean(FOV_NADW_LR[time_present_start:time_present_end])
FOV_ABW_mean_LR		= np.mean(FOV_ABW_LR[time_present_start:time_present_end])
AMOC_mean_LR		= np.mean(AMOC_LR[time_present_start:time_present_end])

#Determine the trends (per century)
FOV_trend_HR, b		= np.polyfit(time_HR[time_future_start:time_future_end], FOV_HR[time_future_start:time_future_end], 1) * 100.0
FOV_ASW_trend_HR, b	= np.polyfit(time_HR[time_future_start:time_future_end], FOV_ASW_HR[time_future_start:time_future_end], 1) * 100.0
FOV_AIW_trend_HR, b	= np.polyfit(time_HR[time_future_start:time_future_end], FOV_AIW_HR[time_future_start:time_future_end], 1) * 100.0
FOV_NADW_trend_HR, b	= np.polyfit(time_HR[time_future_start:time_future_end], FOV_NADW_HR[time_future_start:time_future_end], 1) * 100.0
FOV_ABW_trend_HR, b	= np.polyfit(time_HR[time_future_start:time_future_end], FOV_ABW_HR[time_future_start:time_future_end], 1) * 100.0
AMOC_trend_HR, b	= np.polyfit(time_HR[time_future_start:time_future_end], AMOC_HR[time_future_start:time_future_end], 1) * 100.0

FOV_trend_LR, b		= np.polyfit(time_LR[time_future_start:time_future_end], FOV_LR[time_future_start:time_future_end], 1) * 100.0
FOV_ASW_trend_LR, b	= np.polyfit(time_LR[time_future_start:time_future_end], FOV_ASW_LR[time_future_start:time_future_end], 1) * 100.0
FOV_AIW_trend_LR, b	= np.polyfit(time_LR[time_future_start:time_future_end], FOV_AIW_LR[time_future_start:time_future_end], 1) * 100.0
FOV_NADW_trend_LR, b	= np.polyfit(time_LR[time_future_start:time_future_end], FOV_NADW_LR[time_future_start:time_future_end], 1) * 100.0
FOV_ABW_trend_LR, b	= np.polyfit(time_LR[time_future_start:time_future_end], FOV_ABW_LR[time_future_start:time_future_end], 1) * 100.0
AMOC_trend_LR, b	= np.polyfit(time_LR[time_future_start:time_future_end], AMOC_LR[time_future_start:time_future_end], 1) * 100.0

#-----------------------------------------------------------------------------------------
#Read in the time series of the reanalysis product (1994 - 2020)
time_rean, FOV_rean, FOV_ASW_rean, FOV_AIW_rean, FOV_NADW_rean, FOV_ABW_rean = ReadinData(directory_reanalysis+'FOV_index_section_34S.nc')
time_rean, AMOC_rean							     = ReadinDataAMOC(directory_reanalysis+'AMOC_transport_depth_'+str(depth_min)+'-'+str(depth_max)+'m.nc')

#Determine the present-day time mean
FOV_mean_rean		= np.mean(FOV_rean)
FOV_ASW_mean_rean	= np.mean(FOV_ASW_rean)
FOV_AIW_mean_rean	= np.mean(FOV_AIW_rean)
FOV_NADW_mean_rean	= np.mean(FOV_NADW_rean)
FOV_ABW_mean_rean	= np.mean(FOV_ABW_rean)
AMOC_mean_rean		= np.mean(AMOC_rean)

#-----------------------------------------------------------------------------------------
#Now read in the root-mean-square error w.r.t. reanalysis
fh = netcdf.Dataset(directory_CMIP6+'SALT_VVEL_RMS.nc', 'r')

vel_RMS_all	    = fh.variables['VVEL_RMS'][:] 
salt_RMS_all	= fh.variables['SALT_RMS'][:] 

fh.close()

#-----------------------------------------------------------------------------------------
fh = netcdf.Dataset(directory_iHESP+'SALT_VVEL_RMS_RCP_HR.nc', 'r')

vel_RMS_HR	= fh.variables['VVEL_RMS'][0] 
salt_RMS_HR	= fh.variables['SALT_RMS'][0] 

fh.close()

#-----------------------------------------------------------------------------------------

fh = netcdf.Dataset(directory_iHESP+'SALT_VVEL_RMS_RCP_LR.nc', 'r')

vel_RMS_LR	= fh.variables['VVEL_RMS'][0] 
salt_RMS_LR	= fh.variables['SALT_RMS'][0] 

fh.close()

#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 0], AMOC_all[:, 0], 1)
x_plot	= np.array([transport_salt_all[:, 0].min(), transport_salt_all[:, 0].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((AMOC_all[:, 0] - (a * transport_salt_all[:, 0] + b))**2.0) / np.sum((AMOC_all[:, 0] - np.mean(AMOC_all[:, 0]))**2.0)

print('CMIP6 model regression (FOV vs. AMOC):', a)

#-----------------------------------------------------------------------------------------
fig, ax	= subplots()

ax.fill_between([-0.28, -0.05], y1 = -100, y2 = 100, alpha=0.25, edgecolor='orange', facecolor='orange')
ax.fill_between([-1, 1], y1 = 16, y2 = 19, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 0], AMOC_all[models_1, 0], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 0], AMOC_all[models_2, 0], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 0], AMOC_all[models_3, 0], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.set_xlim(-0.51, 0.51)
ax.set_ylim(0, 30)
ax.set_xlabel('Freshwater transport (Sv)')
ax.set_ylabel('AMOC strength (Sv)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 0]), np.mean(AMOC_all[:, 0]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_mean_HR, AMOC_mean_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_mean_LR, AMOC_mean_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)
ax.scatter(FOV_mean_rean, AMOC_mean_rean, s = 70, color = 'c', label = 'Reanalysis', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 0]), np.mean(AMOC_all[:, 0]), xerr = np.std(transport_salt_all[:, 0]), yerr = np.std(AMOC_all[:, 0]), color = 'k') 
ax.errorbar(FOV_mean_rean, AMOC_mean_rean, xerr = np.std(FOV_rean), yerr = np.std(AMOC_rean), color = 'c') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('a) $F_{\mathrm{OV}}$ (34$^{\circ}$S) and AMOC strength (26$^{\circ}$N), 1994 - 2020')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.53, 0.17, 0.34, 0.2])

ax2.plot(time, np.mean(transport_salt_all_time, axis = 1), '-k')
ax2.fill_between(time, y1 = np.percentile(transport_salt_all_time, 2.5, axis = 1), y2 = np.percentile(transport_salt_all_time, 97.5, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.fill_between(time, y1 = np.percentile(transport_salt_all_time, 25, axis = 1), y2 = np.percentile(transport_salt_all_time, 75, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.set_xlim(1994, 2100)
ax2.set_ylim(-0.5, 0.5)
ax2.grid()
ax2.set_title('$F_{\mathrm{OV}}$ (Sv) in CMIP6', fontsize = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 0], transport_salt_ASW_all[:, 0], 1)
x_plot	= np.array([transport_salt_all[:, 0].min(), transport_salt_all[:, 0].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_ASW_all[:, 0] - (a * transport_salt_all[:, 0] + b))**2.0) / np.sum((transport_salt_ASW_all[:, 0] - np.mean(transport_salt_ASW_all[:, 0]))**2.0)

print('CMIP6 model regression (FOV vs. ASW):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-0.28, -0.05], y1 = -1, y2 = 1, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 0], transport_salt_ASW_all[models_1, 0], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 0], transport_salt_ASW_all[models_2, 0], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 0], transport_salt_ASW_all[models_3, 0], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.51, 0.51)
ax.set_xlabel('Freshwater transport (Sv)')
ax.set_ylabel('ASW freshwater transport (Sv)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_ASW_all[:, 0]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_mean_HR, FOV_ASW_mean_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_mean_LR, FOV_ASW_mean_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)
ax.scatter(FOV_mean_rean, FOV_ASW_mean_rean, s = 70, color = 'c', label = 'Reanalysis', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_ASW_all[:, 0]), xerr = np.std(transport_salt_all[:, 0]), yerr = np.std(transport_salt_ASW_all[:, 0]), color = 'k') 
ax.errorbar(FOV_mean_rean, FOV_ASW_mean_rean, xerr = np.std(FOV_rean), yerr = np.std(FOV_ASW_rean), color = 'c') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('c) $F_{\mathrm{OV}}$ and ASW at 34$^{\circ}$S, 1994 - 2020')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.53, 0.17, 0.34, 0.2])

ax2.plot(time, np.mean(transport_salt_ASW_time, axis = 1), '-k')
ax2.fill_between(time, y1 = np.percentile(transport_salt_ASW_time, 2.5, axis = 1), y2 = np.percentile(transport_salt_ASW_time, 97.5, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.fill_between(time, y1 = np.percentile(transport_salt_ASW_time, 25, axis = 1), y2 = np.percentile(transport_salt_ASW_time, 75, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.set_xlim(1994, 2100)
ax2.set_ylim(-0.3, 0.3)
ax2.set_yticks([-0.3, 0, 0.3])
ax2.grid()
ax2.set_title('ASW (Sv) in CMIP6', fontsize = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 0], transport_salt_AIW_all[:, 0], 1)
x_plot	= np.array([transport_salt_all[:, 0].min(), transport_salt_all[:, 0].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_AIW_all[:, 0] - (a * transport_salt_all[:, 0] + b))**2.0) / np.sum((transport_salt_AIW_all[:, 0] - np.mean(transport_salt_AIW_all[:, 0]))**2.0)

print('CMIP6 model regression (FOV vs. AIW):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-0.28, -0.05], y1 = -1, y2 = 1, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 0], transport_salt_AIW_all[models_1, 0], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 0], transport_salt_AIW_all[models_2, 0], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 0], transport_salt_AIW_all[models_3, 0], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.51, 0.51)
ax.set_xlabel('Freshwater transport (Sv)')
ax.set_ylabel('AIW freshwater transport (Sv)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_AIW_all[:, 0]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_mean_HR, FOV_AIW_mean_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_mean_LR, FOV_AIW_mean_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)
ax.scatter(FOV_mean_rean, FOV_AIW_mean_rean, s = 70, color = 'c', label = 'Reanalysis', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_AIW_all[:, 0]), xerr = np.std(transport_salt_all[:, 0]), yerr = np.std(transport_salt_AIW_all[:, 0]), color = 'k') 
ax.errorbar(FOV_mean_rean, FOV_AIW_mean_rean, xerr = np.std(FOV_rean), yerr = np.std(FOV_AIW_rean), color = 'c') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('e) $F_{\mathrm{OV}}$ and AIW at 34$^{\circ}$S, 1994 - 2020')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.53, 0.17, 0.34, 0.2])

ax2.plot(time, np.mean(transport_salt_AIW_time, axis = 1), '-k')
ax2.fill_between(time, y1 = np.percentile(transport_salt_AIW_time, 2.5, axis = 1), y2 = np.percentile(transport_salt_AIW_time, 97.5, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.fill_between(time, y1 = np.percentile(transport_salt_AIW_time, 25, axis = 1), y2 = np.percentile(transport_salt_AIW_time, 75, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.set_xlim(1994, 2100)
ax2.set_ylim(-0.2, 0.2)
ax2.set_yticks([-0.2, 0, 0.2])
ax2.grid()
ax2.set_title('AIW (Sv) in CMIP6', fontsize = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 0], transport_salt_NADW_all[:, 0], 1)
x_plot	= np.array([transport_salt_all[:, 0].min(), transport_salt_all[:, 0].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_NADW_all[:, 0] - (a * transport_salt_all[:, 0] + b))**2.0) / np.sum((transport_salt_NADW_all[:, 0] - np.mean(transport_salt_NADW_all[:, 0]))**2.0)

print('CMIP6 model regression (FOV vs. NADW):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.fill_between([-0.28, -0.05], y1 = -1, y2 = 1, alpha=0.25, edgecolor='orange', facecolor='orange')

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 0], transport_salt_NADW_all[models_1, 0], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 0], transport_salt_NADW_all[models_2, 0], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 0], transport_salt_NADW_all[models_3, 0], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.51, 0.51)
ax.set_xlabel('Freshwater transport (Sv)')
ax.set_ylabel('NADW freshwater transport (Sv)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_NADW_all[:, 0]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_mean_HR, FOV_NADW_mean_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_mean_LR, FOV_NADW_mean_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)
ax.scatter(FOV_mean_rean, FOV_NADW_mean_rean, s = 70, color = 'c', label = 'Reanalysis', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 0]), np.mean(transport_salt_NADW_all[:, 0]), xerr = np.std(transport_salt_all[:, 0]), yerr = np.std(transport_salt_NADW_all[:, 0]), color = 'k') 
ax.errorbar(FOV_mean_rean, FOV_NADW_mean_rean, xerr = np.std(FOV_rean), yerr = np.std(FOV_NADW_rean), color = 'c') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('g) $F_{\mathrm{OV}}$ and NADW at 34$^{\circ}$S, 1994 - 2020')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.53, 0.17, 0.34, 0.2])

ax2.plot(time, np.mean(transport_salt_NADW_time, axis = 1), '-k')
ax2.fill_between(time, y1 = np.percentile(transport_salt_NADW_time, 2.5, axis = 1), y2 = np.percentile(transport_salt_NADW_time, 97.5, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.fill_between(time, y1 = np.percentile(transport_salt_NADW_time, 25, axis = 1), y2 = np.percentile(transport_salt_NADW_time, 75, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.set_xlim(1994, 2100)
ax2.set_ylim(-0.2, 0.2)
ax2.grid()
ax2.set_title('NADW (Sv) in CMIP6', fontsize = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 1], AMOC_all[:, 1], 1)
x_plot	= np.array([transport_salt_all[:, 1].min(), transport_salt_all[:, 1].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((AMOC_all[:, 1] - (a * transport_salt_all[:, 1] + b))**2.0) / np.sum((AMOC_all[:, 1] - np.mean(AMOC_all[:, 1]))**2.0)

print('CMIP6 model regression (FOV vs. AMOC trends):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 1], AMOC_all[models_1, 1], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 1], AMOC_all[models_2, 1], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 1], AMOC_all[models_3, 1], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)

ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-20, 0)
ax.set_xlabel('Freshwater transport trend (Sv per century)')
ax.set_ylabel('AMOC strength trend (Sv per century)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 1]), np.mean(AMOC_all[:, 1]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_trend_HR, AMOC_trend_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_trend_LR, AMOC_trend_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 1]), np.mean(AMOC_all[:, 1]), xerr = np.std(transport_salt_all[:, 1]), yerr = np.std(AMOC_all[:, 1]), color = 'k') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('b) $F_{\mathrm{OV}}$ and AMOC trends, 2000 - 2100')

#-----------------------------------------------------------------------------------------

ax2 	= fig.add_axes([0.53, 0.17, 0.34, 0.2])

ax2.plot(time, np.mean(AMOC_all_time, axis = 1), '-k')
ax2.fill_between(time, y1 = np.percentile(AMOC_all_time, 2.5, axis = 1), y2 = np.percentile(AMOC_all_time, 97.5, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.fill_between(time, y1 = np.percentile(AMOC_all_time, 25, axis = 1), y2 = np.percentile(AMOC_all_time, 75, axis = 1), alpha=0.25, edgecolor='black', facecolor='black')
ax2.set_xlim(1994, 2100)
ax2.set_ylim(0, 30)
ax2.grid()
ax2.set_title('AMOC strength (Sv) in CMIP6', fontsize = 10)

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 1], transport_salt_ASW_all[:, 1], 1)
x_plot	= np.array([transport_salt_all[:, 1].min(), transport_salt_all[:, 1].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_ASW_all[:, 1] - (a * transport_salt_all[:, 1] + b))**2.0) / np.sum((transport_salt_ASW_all[:, 1] - np.mean(transport_salt_ASW_all[:, 1]))**2.0)

print('CMIP6 model regression (FOV vs. ASW trends):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 1], transport_salt_ASW_all[models_1, 1], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 1], transport_salt_ASW_all[models_2, 1], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 1], transport_salt_ASW_all[models_3, 1], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)

ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.3, 0.3)
ax.set_xlabel('Freshwater transport trend (Sv per century)')
ax.set_ylabel('ASW freshwater transport trend (Sv per century)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_ASW_all[:, 1]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_trend_HR, FOV_ASW_trend_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_trend_LR, FOV_ASW_trend_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_ASW_all[:, 1]), xerr = np.std(transport_salt_all[:, 1]), yerr = np.std(transport_salt_ASW_all[:, 1]), color = 'k') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('d) $F_{\mathrm{OV}}$ and ASW trends, 2000 - 2100')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 1], transport_salt_AIW_all[:, 1], 1)
x_plot	= np.array([transport_salt_all[:, 1].min(), transport_salt_all[:, 1].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_AIW_all[:, 1] - (a * transport_salt_all[:, 1] + b))**2.0) / np.sum((transport_salt_AIW_all[:, 1] - np.mean(transport_salt_AIW_all[:, 1]))**2.0)

print('CMIP6 model regression (FOV vs. AIW trends):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 1], transport_salt_AIW_all[models_1, 1], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 1], transport_salt_AIW_all[models_2, 1], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 1], transport_salt_AIW_all[models_3, 1], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)

ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.3, 0.3)
ax.set_xlabel('Freshwater transport trend (Sv per century)')
ax.set_ylabel('AIW freshwater transport trend (Sv per century)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_AIW_all[:, 1]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_trend_HR, FOV_AIW_trend_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_trend_LR, FOV_AIW_trend_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_AIW_all[:, 1]), xerr = np.std(transport_salt_all[:, 1]), yerr = np.std(transport_salt_AIW_all[:, 1]), color = 'k') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('f) $F_{\mathrm{OV}}$ and AIW trends, 2000 - 2100')

#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------

#Determine the CMIP6 model regression
a, b	= np.polyfit(transport_salt_all[:, 1], transport_salt_NADW_all[:, 1], 1)
x_plot	= np.array([transport_salt_all[:, 1].min(), transport_salt_all[:, 1].max()])
y_plot	= a * x_plot + b
r_2	= 1.0 - np.sum((transport_salt_NADW_all[:, 1] - (a * transport_salt_all[:, 1] + b))**2.0) / np.sum((transport_salt_NADW_all[:, 1] - np.mean(transport_salt_NADW_all[:, 1]))**2.0)

print('CMIP6 model regression (FOV vs. NADW trends):', a)

#-----------------------------------------------------------------------------------------

fig, ax	= subplots()

ax.plot(x_plot, y_plot, '--k', zorder = 0)
ax.scatter(transport_salt_all[models_1, 1], transport_salt_NADW_all[models_1, 1], marker = 'D', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_2, 1], transport_salt_NADW_all[models_2, 1], marker = 'o', color = 'k', s = 50, alpha = 0.5, zorder = 1)
ax.scatter(transport_salt_all[models_3, 1], transport_salt_NADW_all[models_3, 1], marker = 'x', color = 'k', s = 50, alpha = 0.5, zorder = 1)

ax.set_xlim(-0.51, 0.51)
ax.set_ylim(-0.3, 0.3)
ax.set_xlabel('Freshwater transport trend (Sv per century)')
ax.set_ylabel('NADW freshwater transport trend (Sv per century)')
ax.grid()

ax.scatter(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_NADW_all[:, 1]), s = 70, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax.scatter(FOV_trend_HR, FOV_NADW_trend_HR, s = 70, color = 'r', label = 'HR-CESM', zorder = 10)
ax.scatter(FOV_trend_LR, FOV_NADW_trend_LR, s = 70, color = 'b', label = 'LR-CESM', zorder = 10)

ax.errorbar(np.mean(transport_salt_all[:, 1]), np.mean(transport_salt_NADW_all[:, 1]), xerr = np.std(transport_salt_all[:, 1]), yerr = np.std(transport_salt_NADW_all[:, 1]), color = 'k') 

ax.text(0.97, 0.95, '$R^2$ = '+str(round(r_2, 2)), verticalalignment='center', horizontalalignment='right', color = 'k', fontsize = 12, transform = ax.transAxes)

ax.legend(loc='upper left', fancybox=True, shadow=False, scatterpoints=1, ncol = 1, framealpha = 1.0)
ax.set_title('h) $F_{\mathrm{OV}}$ and NADW trends, 2000 - 2100')

#-----------------------------------------------------------------------------------------
ax2 	= fig.add_axes([0.53, 0.2, 0.44, 0.25])

ax2.scatter(salt_RMS_all[models_1], vel_RMS_all[models_1], marker = 'D', color = 'k', s = 20, alpha = 0.5, zorder = 1)
ax2.scatter(salt_RMS_all[models_2], vel_RMS_all[models_2], marker = 'o', color = 'k', s = 20, alpha = 0.5, zorder = 1)
ax2.scatter(salt_RMS_all[models_3], vel_RMS_all[models_3], marker = 'x', color = 'k', s =20, alpha = 0.5, zorder = 1)
ax2.set_xlim(-0.02, 0.6)
ax2.set_ylim(-0.01, 0.15)
ax2.set_xticks([0, 0.2, 0.4, 0.6])
ax2.set_xlabel('Salinity RMS (g kg$^{-1}$)', fontsize = 7.5)
ax2.set_ylabel('Velocity RMS (cm s$^{-1}$)', fontsize = 7.5)
ax2.grid()

ax2.scatter(np.mean(salt_RMS_all), np.mean(vel_RMS_all), s = 30, color = 'k', label = 'CMIP6 mean', zorder = 10)
ax2.scatter(salt_RMS_HR, vel_RMS_HR, s = 30, color = 'r', label = 'HR-CESM', zorder = 10)
ax2.scatter(salt_RMS_LR, vel_RMS_LR, s = 30, color = 'b', label = 'LR-CESM', zorder = 10)
ax2.scatter(0, 0, s = 30, color = 'c', label = 'Reanalysis', zorder = 10)

ax2.errorbar(np.mean(salt_RMS_all), np.mean(vel_RMS_all), xerr = np.std(salt_RMS_all), yerr = np.std(vel_RMS_all), color = 'k') 
ax2.set_title('Salinity and meridional velocity RMS', fontsize = 10)

show()