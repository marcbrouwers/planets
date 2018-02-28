# -*- coding: utf-8 -*-
#file pars.py

import math as m

#constants
G = 6.67408e-08 #Newton’s gravitational constant
mSun= 1.9884754153381438e+33 #mass of the Sun (in grams)
mEarth = 5.9722e+27
au = 1.495978707e13 #astronomical unit (in cm)
yr = 2*m.pi /m.sqrt(G*mSun/au**3) #1 year in seconds
vKep_au = m.sqrt(G*mSun /au)


#problem parameters
Np = 1 #2 particles
e_max = 0.
R_dust = 1.