# -*- coding: utf-8 -*-

#main program file
import numpy as np
import functions as fn
import pars
import matplotlib.pyplot as plt


#declarations
tfinal =  1 * pars.yr #end time (seconds)


for dt_yr in [1e-3]:
    xdust, vdust, mdust = fn.init_2body()
    r0 = np.linalg.norm(xdust, axis=1)
    v0 = np.linalg.norm(vdust, axis=1)
    etot0 = fn.energies (xdust, vdust, mdust)
    dt = dt_yr * pars.yr
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6))
    i = 0
    time = 0
    while time<tfinal:
        xdust, vdust = fn.integrator(xdust, vdust, mdust, dt)
        etot = fn.energies (xdust, vdust, mdust)
        time += dt
        
        i+= 1
        if i % 5 == 0:
            #ax1.scatter(time*np.ones(pars.Np), (np.linalg.norm(xdust, axis=1) - r0)/r0, color='k')
            #ax1.scatter(xdust[:,0], xdust[:,1], color='k')
            vr = np.inner(xdust, vdust) / r0
            vt = np.sqrt(np.linalg.norm(vdust, axis=1)**2 - vr**2)
            ax1.scatter(time / (2 * m.pi * pars.au / pars.vKep_au) * np.ones(pars.Np), vt, color='k')
            ax2.scatter(time*np.ones(pars.Np), vr, color='k')
    plt.show()
    etot = fn.energies (xdust, vdust, mdust)
    error = (etot -etot0) /etot0
    #print error