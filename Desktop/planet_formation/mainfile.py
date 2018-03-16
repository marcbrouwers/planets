# -*- coding: utf-8 -*-

#main program file
from __future__ import division
import numpy as np
import functions as fn
import pars
import matplotlib.pyplot as plt
import math as m


#declarations
tfinal =  10 * pars.yr #end time (seconds)


for dt_yr in [1e-3]:
    plot_steps = 200
    plot_interval = int(tfinal / pars.yr / dt_yr / plot_steps)
    xdust, vdust, mdust = fn.init_2body()
#    etot0 = fn.energies (xdust, vdust, mdust)
    dt = dt_yr * pars.yr
    
    fig1, axes = plt.subplots(2,2, figsize = (12,8))
    i = 0
    time = 0
    n = 0 # number of eddies
    while time<tfinal:
        xdust, vdust = fn.integrator(xdust, vdust, mdust, dt, pars.use_eddie)
        
        if time < (n+1)*pars.t_eddie + n * pars.t_no_eddie and pars.use_eddie is True:
            Ed = True
        elif time < (n+1)*pars.t_eddie + (n+1)*pars.t_no_eddie:
            Ed = False
        else:
            n += 1
        time += dt
        
        i+= 1
        if i % plot_interval == 0:
            
            #ax1.scatter(xdust[0,:] * np.cos(xdust[1,:])/pars.au, xdust[0,:] * np.sin(xdust[1,:])/pars.au, color='k')
            
            axes[0,0].scatter(time/pars.tKep_au * np.ones(pars.Np), xdust[0,:], color='k') # vr
            axes[0,1].scatter(time/pars.tKep_au * np.ones(pars.Np), xdust[2,:], color='k') # vz
            axes[1,0].scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[0,:], color='k') # r
            axes[1,1].scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[2,:], color='k') # z
            
            axes[0,0].set_xlabel('t (orbital times)')
            axes[1,0].set_xlabel('t (orbital times)')
            axes[0,1].set_xlabel('t (orbital times)')
            axes[1,1].set_xlabel('t (orbital times)')
            axes[0,0].set_ylabel('r (cm)')
            axes[1,0].set_ylabel('vr (cm/s)')
            axes[0,1].set_ylabel('z (cm)')
            axes[1,1].set_ylabel('vz (cm/s)')
            
    print n, 'total eddies'   
    plt.show()