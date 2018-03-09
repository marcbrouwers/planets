# -*- coding: utf-8 -*-


from __future__ import division
import numpy as np
import pars
import math as m
from timeit import default_timer as timer


def init_2body ():
    '''docstring'''
    
    r = np.ones(pars.Np) * pars.au # starting radii particles
    phi = np.zeros(pars.Np) # starting angles in xy plane
    z = np.zeros(pars.Np) # starting z positions
    xdust = np.stack((r, phi, z), axis=-1) # position vector
    
    T_gas, rho_gas, omega_gas, eta = disk_properties(pars.r0)
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * pars.visc)
    v_r = - eta / (v_kep(r) * t_stop / r +  r / (v_kep(r) * t_stop)) * v_kep(r)
    v_phi = np.sqrt(v_kep(r)**2 + r *v_r / t_stop)
    vdust = np.stack((v_r * np.ones(pars.Np), v_phi/r, np.zeros(pars.Np)), axis=-1) # position vector
    mdust = np.ones(pars.Np)
    
    return xdust, vdust, mdust


def v_kep (r):
    return np.sqrt(pars.G * pars.mSun / r)


def disk_properties (r):
    '''docstring'''
    
    T_disk = 150 * (r/(5.2*pars.au))**-.5
    rho_disk = 5e-11 * (r/(5.2*pars.au))**-2.75
    cs_disk = pars.cs_cst * np.sqrt(T_disk)
    eta = pars.n * cs_disk**2 /  v_kep(r)**2
    
    rho_gas = 0.01*rho_disk
    omega_gas = v_kep(r) / r * np.sqrt(1. - eta)
    T_gas = T_disk
    
    return T_gas, rho_gas, omega_gas, eta
    

def accelerations (xobj, vobj, mobj):
    '''docstring'''
    
    T_gas, rho_gas, omega_gas, eta = disk_properties(xobj[:,0])
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * pars.visc)
    
    vgas = np.stack((np.zeros(pars.Np), np.ones(pars.Np) * omega_gas, np.zeros(pars.Np)), axis=-1)
    
    acc_r_Gravity = - pars.G * pars.mSun / xobj[:,0]**2
    acc_r_centrifical = vobj[:,1]**2 * xobj[:,0]
    acc_phi_centrifical = - 2*vobj[:,0] / xobj[:,0] * vobj[:,1]
    acc_drag = (vgas - vobj) / t_stop
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
