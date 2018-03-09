# -*- coding: utf-8 -*-



import numpy as np
import pars
import math as m
from timeit import default_timer as timer


def init_2body ():
    '''docstring'''
    
    r = np.ones(pars.Np) * pars.au # starting radii particles
    phi = np.zeros(pars.Np) # starting angles in xy plane
    z = np.zeros(pars.Np) # starting z positions
    v = np.sqrt(pars.G*pars.mSun/r) # keplerian starting velocity
    xdust = np.stack((r, phi, z), axis=-1) # position vector
    vdust = 0.9999*np.stack((np.zeros(pars.Np),v/r, np.zeros(pars.Np)), axis=-1) # velocity vector
    mdust = np.ones(pars.Np)
    
    return xdust, vdust, mdust


def v_kep (xobj):
    return np.sqrt(pars.G * pars.mSun / xobj[:,0])
    
def disk_properties (r):
    '''docstring'''
    
    T_disk = 150 * (r/(5.2*pars.au))**-.5
    rho_disk = 5e-11 * (r/(5.2*pars.au))**-2.75
    
    return T_disk, rho_disk
    

def accelerations (xobj, vobj, mobj):
    '''docstring'''
    
    tstop = 1 * pars.yr
    vgas = np.stack((np.zeros(pars.Np), 0.996*vobj[:,1], np.zeros(pars.Np)), axis=-1)
    acc_r_Gravity = - pars.G * pars.mSun / xobj[:,0]**2
    acc_r_centrifical = vobj[:,1]**2 * xobj[:,0]
    acc_phi_centrifical = - 2*vobj[:,0] / xobj[:,0] * vobj[:,1]
    acc_drag = (vgas - vobj) / tstop
    acc_total = acc_drag
    acc_total[:,0] += acc_r_Gravity + acc_r_centrifical
    acc_total[:,1] += acc_phi_centrifical
    
    return acc_total


# =============================================================================
# def energies (xobj, vobj, mobj):
#     '''docstring'''
#     
#     E_kin = 0.5*mobj*np.linalg.norm(vobj,axis=1)**2
#     E_pot = - (pars.G*pars.mSun*mobj/np.linalg.norm(xobj,axis=1))
# 
#     E_tot = E_kin[0] + E_pot[0]
# 
#     return E_tot
# =============================================================================


def integrator (xarr, varr, marr, dt):
    '''docstring'''
    
    acc = accelerations(xarr, varr, marr)
    varr += acc*dt/2
    xarr += varr*dt
    accnew = accelerations(xarr, varr, marr) #update accelerations
    varr += accnew*dt/2

# =============================================================================
#     a = accelerations(xarr,varr, marr)
#     v2 = varr + a*dt*.5
#     a2 = accelerations(xarr + varr*dt*.5, varr +a*dt*.5, marr)
#     v3 = varr + a2 * dt*.5
#     a3 = accelerations(xarr + v2*dt*.5, varr +a*dt*.5, marr)
#     v4 = varr + a3 * dt
#     a4 = accelerations(xarr + v3*dt, varr +a*dt*.5, marr)
#     xarr += (varr + 2*v2 + 2*v3 + v4)*dt / 6
#     varr += (a + 2*a2 + 2*a3 + a4)*dt / 6
# =============================================================================

    return xarr, varr
