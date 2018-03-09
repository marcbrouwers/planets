# -*- coding: utf-8 -*-
#file pars.py

from __future__ import division
import math as m

#constants
G = 6.67408e-08 #Newton’s gravitational constant
mSun= 1.9884754153381438e+33 #mass of the Sun (in grams)
mEarth = 5.9722e+27
au = 1.495978707e13 #astronomical unit (in cm)
yr = 2*m.pi /m.sqrt(G*mSun/au**3) #1 year in seconds
vKep_au = m.sqrt(G*mSun /au)
tKep_au = 2 * m.pi * au / vKep_au

X, Y = 0.75, 0.25
u =  1.660539040e-24
m_H2 = 2.*u
m_He = 4.*u
k = 1.3806503e-16
gamma_av = 1.4
m_av = 1. / (X / m_H2 + Y / m_He)
cs_cst = m.sqrt(gamma_av * k / m_av)
n = 3

#problem parameters
Np = 1 #2 particles
e_max = 0.
R_dust = 100.
rho_dust = 3.
visc = 2.4e+5
r0 = 1 * au