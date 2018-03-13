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
    z = r * 0.01 # starting z positions
    xdust = np.stack((r, phi, z), axis=0) # position vector
    
    T_gas, rho_gas, omega_gas, eta = disk_properties(r, z)
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * pars.visc)
    print t_stop / pars.yr
    v_r = - eta / (v_kep(r) * t_stop / r +  r / (v_kep(r) * t_stop)) * v_kep(r)
    v_phi = np.sqrt(v_kep(r)**2 + r *v_r / t_stop)
    v_z = -omega_gas**2 * t_stop * z
    print v_z
    vdust = np.stack((v_r * np.ones(pars.Np), v_phi/r, v_z), axis=0) # position vector
    mdust = np.ones(pars.Np)
    
    return xdust, vdust, mdust

def v_kep (r):
    return np.sqrt(pars.G * pars.mSun / r)


def disk_properties (r, z):
    '''docstring'''
    
    T = 600 * (r/pars.au)**-.5
    sigma = 1.7e+3 * (r/pars.au)**-1.5
    cs = pars.cs_cst * np.sqrt(T)
    omega = v_kep(r) / r
    h = cs / omega
    eta = pars.n * cs**2 /  v_kep(r)**2
    rho_disk = sigma / (np.sqrt(2*m.pi)*h)
    
    rho_gas = 0.01*rho_disk*np.exp(-0.5*(z/h)**2)
    omega_gas = omega * np.sqrt(1. - eta)

    return T, rho_gas, omega_gas, eta
    

def accelerations (xobj, vobj, mobj):
    '''docstring'''
    
    r, phi, z = xobj[0,:], xobj[1,:], xobj[2,:]
    v_r, omega_phi, v_z = vobj[0,:], vobj[1,:], vobj[2,:]
    theta = np.arctan(z/r)
    T_gas, rho_gas, omega_gas, eta = disk_properties(r, z)
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * pars.visc)
    
    vgas = np.stack((np.zeros(pars.Np), np.ones(pars.Np) * omega_gas, np.zeros(pars.Np)), axis=0)
    
    acc_Gravity_star = - pars.G * pars.mSun / (r**2 + z**2)
    acc_r_centrifical = omega_phi**2 * r
    acc_phi_centrifical = - 2*v_r / r * omega_phi
    acc_drag = (vgas - vobj) / t_stop
    acc_total = acc_drag
    acc_total[0,:] += acc_Gravity_star * np.cos(theta)+ acc_r_centrifical
    acc_total[1,:] += acc_phi_centrifical
    acc_total[2,:] += acc_Gravity_star * np.sin(theta)

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
