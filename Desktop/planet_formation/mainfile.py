# -*- coding: utf-8 -*-

#main program file
import numpy as np
import functions as fn
import pars
import matplotlib.pyplot as plt
import math as m


#declarations
tfinal =  1.7 * pars.yr #end time (seconds)


for dt_yr in [1e-3]:
    xdust, vdust, mdust = fn.init_2body()
#    etot0 = fn.energies (xdust, vdust, mdust)
    dt = dt_yr * pars.yr
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6))
    i = 0
    time = 0
    while time<tfinal:
        xdust, vdust = fn.integrator(xdust, vdust, mdust, dt)
#        etot = fn.energies (xdust, vdust, mdust)
        time += dt
        
        i+= 1
        if i % 10 == 0:
            ax1.scatter(time/pars.tKep_au * np.ones(pars.Np), vdust[:, 1] * xdust[:,0], color='k')
            ax2.scatter(xdust[:,0] * np.cos(xdust[:,1])/pars.au, xdust[:,0] * np.sin(xdust[:,1])/pars.au, color='k')
            #ax1.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), vt, color='k')
            #ax2.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), vr, color='k')
            #ax1.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), force_plot, color='k')
    plt.show()
    
#    etot = fn.energies (xdust, vdust, mdust)
#    error = (etot -etot0) /etot0
#    print error