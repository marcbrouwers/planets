# -*- coding: utf-8 -*-



import numpy as np
import pars
import math as m
from timeit import default_timer as timer



def init_2body ():

    """
    construct the 2body problem; initialize at aphelion
    """
    r = np.random.uniform(1, high=1, size=pars.Np) * pars.au
    ecc = np.random.uniform(0., high=pars.e_max, size=pars.Np)
    v = np.sqrt(pars.G*pars.mSun/r)*np.sqrt((1-ecc) /(1+ecc))
    xdust = np.stack((r,np.zeros(pars.Np), np.zeros(pars.Np)), axis=-1)
    vdust = np.stack((np.zeros(pars.Np),v, np.zeros(pars.Np)), axis=-1)    
    mdust = np.ones((pars.Np,1))
    return xdust, vdust, mdust

def disk_properties (r):
    T_disk = 150 * (r/(5.2*pars.au))**-.5
    rho_disk = 5e-11 * (r/(5.2*pars.au))**-2.75
    
    return rho_disk
    

def forces (xobj, vobj, mobj):

    # Gravitational force
    r_obj = np.linalg.norm(xobj,axis=-1)[:,None]
    force_G = - pars.G * pars.mSun * mobj * xobj / r_obj**3
    tstop = 10 * pars.yr
    force_D = - 0.01 * vobj / tstop
    
    force = force_G + force_D
    return force


def energies (xobj, vobj, mobj):
    
    E_kin = 0.5*mobj*np.linalg.norm(vobj,axis=1)**2
    E_pot = - (pars.G*pars.mSun*mobj/np.linalg.norm(xobj,axis=1))

    E_tot = E_kin[0] + E_pot[0]

    return E_tot


def integrator (xarr, varr, marr, dt):

    acc = forces(xarr, varr, marr) / marr
    varr += acc*dt/2
    xarr += varr*dt
    accnew = forces(xarr, varr, marr) / marr #update accelerations
    varr += accnew*dt/2

# =============================================================================
#     a = forces(xarr,varr, marr) / marr
#     v2 = varr + a*dt*.5
#     a2 = forces(xarr + varr*dt*.5, varr +a*dt*.5, marr) / marr
#     v3 = varr + a2 * dt*.5
#     a3 = forces(xarr + v2*dt*.5, varr +a*dt*.5, marr) / marr
#     v4 = varr + a3 * dt
#     a4 = forces(xarr + v3*dt, varr +a*dt*.5, marr) / marr
#     xarr += (varr + 2*v2 + 2*v3 + v4)*dt / 6
#     varr += (a + 2*a2 + 2*a3 + a4)*dt / 6
# =============================================================================

    return xarr, varr
