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
    z = r * 0. # starting z positions
    xdust = np.stack((r, phi, z), axis=0) # position vector
    
    T_gas, rho_gas, omega_gas, eta, v_th = disk_properties(r, z)
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * v_th)
    vkep = v_kep(r, z)
    v_r = - eta / (vkep * t_stop / r +  r / (vkep * t_stop)) * vkep
    v_phi = np.sqrt(vkep**2 + r *v_r / t_stop)
    v_z = -omega_gas**2 * t_stop * z
    vdust = np.stack((v_r, v_phi/r, v_z), axis=0) # position vector
    
    mdust = 4/3*m.pi*pars.R_dust**3
    
    print t_stop/pars.yr
    
    return xdust, vdust, mdust


def v_kep (r, z):
    theta = np.arctan(r/z)
    
    return np.sqrt(pars.G * pars.mSun) * np.sqrt(r / (r**2 + z**2)) * np.sin(theta)


def disk_properties (r, z):
    '''docstring'''
    
    T = 600 * (r/pars.au)**-.5
   
    sigma = 1.7e+3 * (r/pars.au)**-1.5
    cs = pars.cs_cst * np.sqrt(T)
    omega = v_kep(r,z) / r
    h = cs / omega
    eta = pars.n * cs**2 /  v_kep(r,z)**2
    rho_disk = sigma / (2*h)
    
    rho_gas = 0.01*rho_disk*np.exp(-0.5*(z/h)**2)
    omega_gas = omega * np.sqrt(1. - eta)
    v_th = m.sqrt(8/m.pi)*cs
    
    return T, rho_gas, omega_gas, eta, v_th
    

def accelerations (xobj, vobj, mobj, Ed):
    '''docstring'''
    
    r, phi, z = xobj[0,:], xobj[1,:], xobj[2,:]
    v_r, omega_phi, v_z = vobj[0,:], vobj[1,:], vobj[2,:]
    theta = np.arctan(z/r)
    T_gas, rho_gas, omega_gas, eta, v_th = disk_properties(r, z)
    t_stop = pars.rho_dust * pars.R_dust / (rho_gas * v_th)

    if Ed is True:
        angle_eddie = 2*m.pi*np.random.rand(pars.Np)[0]
        v_eddie_phi = pars.v_eddie * np.cos(angle_eddie)
        omega_eddie_phi = v_eddie_phi / r
        v_eddie_r = pars.v_eddie * np.sin(angle_eddie)
    else:
        omega_eddie_phi, v_eddie_r = np.zeros(pars.Np), np.zeros(pars.Np)
    
    vgas = np.stack((np.zeros(pars.Np) + v_eddie_r, np.ones(pars.Np) * omega_gas + omega_eddie_phi, np.zeros(pars.Np)), axis=0)
    
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


def integrator (xarr, varr, marr, dt, Ed):
    '''docstring'''
    
    acc = accelerations(xarr, varr, marr, Ed)
    varr += acc*dt/2
    xarr += varr*dt
    accnew = accelerations(xarr, varr, marr, Ed) #update accelerations
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
