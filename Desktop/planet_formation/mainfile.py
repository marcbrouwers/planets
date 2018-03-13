# -*- coding: utf-8 -*-

#main program file
from __future__ import division
import numpy as np
import functions as fn
import pars
import matplotlib.pyplot as plt
import math as m


#declarations
tfinal =  2 * pars.yr #end time (seconds)


for dt_yr in [1e-3]:
    plot_steps = 200
    plot_interval = int(tfinal / pars.yr / dt_yr / plot_steps)
    xdust, vdust, mdust = fn.init_2body()
#    etot0 = fn.energies (xdust, vdust, mdust)
    dt = dt_yr * pars.yr
    
    fig1, axes = plt.subplots(2,2, figsize = (12,8))
    i = 0
    time = 0
    while time<tfinal:
        xdust, vdust = fn.integrator(xdust, vdust, mdust, dt)
#        etot = fn.energies (xdust, vdust, mdust)
        time += dt
        
        i+= 1
        if i % plot_interval == 0:
            
            #ax.scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[1, :] * xdust[0,:], color='k')
            #ax1.scatter(xdust[0,:] * np.cos(xdust[1,:])/pars.au, xdust[0,:] * np.sin(xdust[1,:])/pars.au, color='k')
            axes[0,0].scatter(time/pars.tKep_au * np.ones(pars.Np), xdust[0,:], color='k') # vr
            axes[0,1].scatter(time/pars.tKep_au * np.ones(pars.Np), xdust[2,:], color='k') # vz
            axes[1,0].scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[0,:], color='k') # r
            axes[1,1].scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[2,:], color='k') # z
            #ax2.scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[:,1] * xdust[:,0] - fn.v_kep(xdust[:,0]), color='k')
            #e = -pars.G * pars.mSun / xdust[0,:] + 0.5 * (vdust[0,:]**2 + (vdust[1,:]*xdust[0,:])**2)
            #ax2.scatter(time/pars.tKep_au * np.ones(pars.Np), xdust[0,:], color='k')
            #ax1.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), vt, color='k')
            #ax2.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), vr, color='k')
            #ax1.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), force_plot, color='k')
            axes[0,0].set_xlabel('t (orbital times)')
            axes[1,0].set_xlabel('t (orbital times)')
            axes[0,1].set_xlabel('t (orbital times)')
            axes[1,1].set_xlabel('t (orbital times)')
            axes[0,0].set_ylabel('r (cm/s)')
            axes[1,0].set_ylabel('vr (cm/s)')
            axes[0,1].set_ylabel('z (cm/s)')
            axes[1,1].set_ylabel('vz (cm/s)')
    plt.show()
    
#    etot = fn.energies (xdust, vdust, mdust)
#    error = (etot -etot0) /etot0
#    print error