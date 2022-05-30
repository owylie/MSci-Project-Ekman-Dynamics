import numpy as np
import datetime as dt
from dateutil.relativedelta import relativedelta
import matplotlib.pyplot as plt
import matplotlib.colors as cl
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import glob
from sys import path
from os.path import exists
from imp import reload
from netCDF4 import Dataset
# path.insert(0, '/Users/H/WAVES/geo_data_group/')
path.insert(0, '/Users/owylie/MSci/')
import grid_set as gs

import Ekman_pumping as ep



#need Pathfinder grid points
m = ccrs.NorthPolarStereo()
GPathfinder = gs.grid_set(m)
GPathfinder.load_grid('/Users/owylie/MSci/Pathfinder_gs.npz')



#m_use = ep.m_AIDJEX_0 #(has angles = 0)
m_use = ep.m_AIDJEX #check Ekman_pumping.py to see what constants this includes
[Ca,Co] = m_use[0:2] #this just gets the values of Ca and Co, can change to regular values if you want

MLD=20.0 #Ekman depth 
omega= 2.*np.pi/24./60./60. #Earths rotation
fcor = 2.*omega*np.sin(np.deg2rad(GPathfinder.lats.T)) #coriolis force




#setting up empty arrays with the size of GPathfinder
ue = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
ueg0= np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
tau_a = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
tau_i = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
tau_i0 = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
tau_g = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
tau_all= np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])

np_ice = 0
np_ocn = 0



### FUNCTION THAT CALCULATES THE TAUS AND THE EKMAN CURRENTS

def current(uwp,vwp,upp,vpp,ugp,vgp,alpha):
    
    
    #setting up empty arrays with the size of GPathfinder
    ue = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    ueg0= np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_a = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_i = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_i0 = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_g = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_all= np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])

    np_ice = 0
    np_ocn = 0
    
    
    for i in range (len(GPathfinder.xpts)):
        for j in range (len(GPathfinder.ypts)):
        
        
            #ICE-COVERED
            if np.isfinite(ugp[i,j]) and np.isfinite(vgp[i,j]) and np.isfinite(alpha[i,j]) and alpha[i,j]> 0.0001:
                #np_ice += 1
            
                #lambda is basically creating a small function, solving x = Ek_u and y = Ek_v
                g = lambda x,y: ep.ueck_solve_quad(m_use,uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],x,y,
                                                   ugp[i,j],vgp[i,j],alpha[i,j],fcor[i,j],MLD)
            
                x0 = g(0, 0) #0,0 are Ekman current u&v vels
            
                #x and y components of minimized pumping force
                ue[i,j,0]=x0[0]
                ue[i,j,1]=x0[1] 
            
                #x and y components of air and ice stress
                [tau_a[i,j,0],tau_a[i,j,1],tau_i[i,j,0],tau_i[i,j,1]] = ep.taus_quad(
                    m_use,uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],ue[i,j,0],ue[i,j,1],ugp[i,j],vgp[i,j],alpha[i,j])
            

                # solve again without geostrophic currents - to see the difference
                g = lambda x,y: ep.ueck_solve_quad(m_use,uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],x,y,
                                                   0.0,0.0,alpha[i,j],fcor[i,j],MLD)
                x0 = g(0, 0)
                ueg0[i,j,0]=x0[0]
                ueg0[i,j,1]=x0[1]
            
                [poo,poo,tau_i0[i,j,0],tau_i0[i,j,1]] = ep.taus_quad(m_use,uwp[i,j],vwp[i,j],
                                                                     upp[i,j],vpp[i,j],ueg0[i,j,0],ueg0[i,j,1],
                                                                     0.0,0.0,alpha[i,j])
            
            
                #x and y components of geostrophic current forcing and total forcing
                tau_g[i,j,0] = tau_i[i,j,0] - tau_i0[i,j,0]
                tau_g[i,j,1] = tau_i[i,j,1] - tau_i0[i,j,1]
                tau_all[i,j,0] = tau_a[i,j,0] + tau_i[i,j,0]
                tau_all[i,j,1] = tau_a[i,j,1] + tau_i[i,j,1] 
            
            #OPEN OCEAN   
            elif np.isfinite(ugp[i,j]) and np.isfinite(vgp[i,j]): #and np.isfinite(uwp[i,j]) and np.isfinite(vwp[i,j]):
                #np_ocn += 1
            
                #x and y components of air stress (= total stress)
                tau_a[i,j,0],tau_a[i,j,1] = ep.taus_quad_open(m_use,uwp[i,j],vwp[i,j])
                tau_all[i,j,:] = tau_a[i,j,:]
            
                ue[i,j,0],ue[i,j,1] = ep.ueck_quad(tau_all[i,j,0],tau_all[i,j,1],fcor[i,j],MLD)
            
                ueg0[i,j,:]=np.nan
                tau_i[i,j,:] = np.nan
                tau_i0[i,j,:] = np.nan
                tau_g[i,j,:] = np.nan

            
            #OVER LAND
            else:
                ue[i,j,:]=np.nan
                ueg0[i,j,:]=np.nan
                tau_all[i,j,:] = np.nan
                tau_a[i,j,:] = np.nan
                tau_i[i,j,:] = np.nan
                tau_i0[i,j,:] = np.nan
                tau_g[i,j,:] = np.nan
        
            #print(tau_all[i,j])
        
    #print("Pump calc complete,",np_ice,"ice points,",np_ocn,"ocn points")
    
    return ue, ueg0, tau_all, tau_a, tau_i, tau_i0, tau_g
    
    
    
    
    

#### THIS IS INCORRECT, DO NOT USE!!!   
def weighted(uwp,vwp,upp,vpp,ugp,vgp,alpha):
    

    
    #setting up empty arrays with the size of GPathfinder
    ue = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    ueg0 = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    #ue_open = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    ue_final = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    tau_a = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    #tau_a_open = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    tau_i = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_i0 = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_g = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    tau_all = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    #tau_ice_covered = np.zeros([len(GPathfinder.xpts),len(GPathfinder.ypts),2])
    #tau_open_ocean = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])
    #tau_total = np.zeros([len(GPathfinder.xpts),len(GPathfinder.xpts),2])


    #np_ice = 0
    #np_ocn = 0
    
    

    for i in range (len(GPathfinder.xpts)):
        for j in range (len(GPathfinder.ypts)):
        
            #OVER OCEAN
            if np.isfinite(ugp[i,j]) and np.isfinite(vgp[i,j]): 
        
        
            
                ########### calculating the ice covered contribution ###########
                if np.isfinite(alpha[i,j]): # and alpha[i,j]>0.00001:
            
                    #lambda is basically creating a small function, solving x = Ek_u and y = Ek_v
                    g = lambda x,y: ep.ueck_solve_quad(m_use,
                                               uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],
                                               x,y,ugp[i,j],vgp[i,j],
                                               alpha[i,j],fcor[i,j],MLD)
            
                    x0 = g(0, 0) #solve for minimum
                    ue[i,j,0]=x0[0] #Ekman vels
                    ue[i,j,1]=x0[1] 
                    [tau_a[i,j,0],tau_a[i,j,1],tau_i[i,j,0],tau_i[i,j,1]] = ep.taus_quad(
                    m_use,uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],ue[i,j,0],ue[i,j,1],
                    ugp[i,j],vgp[i,j],alpha[i,j]) #x and y components of air and ice taus
          
                        ######### (NO GEOS) #########
            
                    g = lambda x,y: ep.ueck_solve_quad(m_use,
                                               uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],
                                               x,y,0.0,0.0,
                                               alpha[i,j],fcor[i,j],MLD)
                    x0 = g(0, 0) #solve for minimum
                    ueg0[i,j,0]=x0[0] #Ekman vels
                    ueg0[i,j,1]=x0[1]
                    [poo,poo,tau_i0[i,j,0],tau_i0[i,j,1]] = ep.taus_quad(
                    m_use,uwp[i,j],vwp[i,j],upp[i,j],vpp[i,j],ueg0[i,j,0],ueg0[i,j,1],
                    0.0,0.0,alpha[i,j]) #x and y components of ice taus (no geo)
            
                    #total taus and geostrophic contribution
                    tau_g[i,j,0] = tau_i[i,j,0] - tau_i0[i,j,0]
                    tau_g[i,j,1] = tau_i[i,j,1] - tau_i0[i,j,1]
                    #tau_ice_covered[i,j,0] = tau_a[i,j,0] + tau_i[i,j,0]
                    #tau_ice_covered[i,j,1] = tau_a[i,j,1] + tau_i[i,j,1] 
                    tau_all[i,j,0] = tau_a[i,j,0] + tau_i[i,j,0]
                    tau_all[i,j,1] = tau_a[i,j,1] + tau_i[i,j,1] 
        
                else:
                    #tau_ice_covered[i,j,:] = 0.0
                    tau_all[i,j,:] = 0.0
                    tau_g[i,j,:] = 0.0
                    ueg0[i,j,:] = 0.0
                    ue[i,j,:] = 0.0
                    tau_i[i,j,:] = 0.0
                    tau_i0[i,j,:] = 0.0
                    tau_a[i,j,:] = 0.0
        
                ############# calculatiing the open ocean contribution ##############
            
                #tau_a_open[i,j,0],tau_a_open[i,j,1] = ep.taus_quad_open(m_use,uwp[i,j],vwp[i,j])
                #tau_open_ocean[i,j,:] = tau_a_open[i,j,:]
            
                #ue_open[i,j,0],ue[i,j,1] = ep.ueck_quad(tau_open_ocean[i,j,0],tau_open_ocean[i,j,1],fcor[i,j],MLD)
        
        
                ######### calculating the weighted average and final currents ##########
            
                #tau_total[i,j,:] = alpha[i,j]*tau_ice_covered[i,j,:] + (1-alpha[i,j])*tau_open_ocean[i,j,:]
                #ue_final[i,j,0],ue_final[i,j,1] = ep.ueck_quad(
                #tau_total[i,j,0],tau_total[i,j,1],fcor[i,j],MLD) #recalculating Ek vels with new tau
            
        
            #OVER LAND
            else:
                ue[i,j,:]=np.nan
                #ue_open[i,j,:]=np.nan
                ueg0[i,j,:]=np.nan
                #ue_final[i,j,:]=np.nan
                tau_a[i,j,:]=np.nan
                #tau_a_open[i,j,:]=np.nan
                tau_i[i,j,:]=np.nan
                tau_i0[i,j,:]=np.nan
                tau_g[i,j,:]=np.nan
                #tau_ice_covered[i,j,:]=np.nan
                #tau_open_ocean[i,j,:]=np.nan
                #tau_total[i,j,:]=np.nan
                
    #return ue_final,tau_total
    return ue, tau_all, tau_a, tau_i, tau_g



    
    
    
    
    
    
### FUNCTION THAT CALCULATES THE PUMPING
def pump(tau_a,tau_i,tau_g,tau_all):
    pump_scale = 60*60*24*365  # m/year
    pump_a = - gs.geo_curl(tau_a[:,:,0],tau_a[:,:,1],GPathfinder)*pump_scale/1027.5/fcor
    pump_i = - gs.geo_curl(tau_i[:,:,0],tau_i[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    pump_g = - gs.geo_curl(tau_g[:,:,0],tau_g[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    pump_all = - gs.geo_curl(tau_all[:,:,0],tau_all[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    
    return pump_a, pump_i, pump_g, pump_all


def weighted_pump(tau_all):
    pump_scale = 60*60*24*365  # m/year
    #pump_a = gs.geo_curl(tau_a[:,:,0],tau_a[:,:,1],GPathfinder)*pump_scale/1027.5/fcor
    #pump_i = gs.geo_curl(tau_i[:,:,0],tau_i[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    #pump_g = gs.geo_curl(tau_g[:,:,0],tau_g[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    pump_all = - gs.geo_curl(tau_all[:,:,0],tau_all[:,:,1],GPathfinder)*pump_scale/ep.rhoo/fcor
    #pump_all2 = pump_a + pump_i
    
    return pump_all