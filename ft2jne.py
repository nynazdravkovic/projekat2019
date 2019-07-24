# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:01:39 2019

@author: nina
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.integrate as integrate
#import sympy
#from sympy import *
#from sympy.solvers.solvers import solve_linear_system

#za resavanje simbolicki 
#dc=sympy.Symbol('dc')
#dp=sympy.Symbol('dp')
#oc=sympy.Symbol('oc')
#op=sympy.Symbol('op')
#ro00=sympy.Symbol('ro00')
#ro01=sympy.Symbol('ro01')
#ro02=sympy.Symbol('ro02')
#gamma01=sympy.Symbol('gamma01')
#gamma02=sympy.Symbol('gamma02')
#w = sympy.Symbol('w')
#I = sympy.I
#za direktno resavanje
dc=0.
dp=0.
oc=0.5
op=0.5
gamma01=1.
gamma02=1.
#w = 1.
I = complex(0,1)
omega = np.linspace(0,100,100)
E = []
#jednacina = sympy.solve([I*w*ro01-I*0.5*oc*ro02+I*0.5*op+ro01*(-gamma01+I*(dc-dp)),I*w*ro02-I*0.5*oc*ro01-I*0.5*op+ro02*(-gamma02 - I*dp)],[ro01,ro02])

def f (w):
    sistem = [[I*w+(-gamma01+I*(dc-dp)),-I*0.5*oc],[-I*0.5*oc,I*w+(-gamma02 - I*dp)]]
    resenje = [-I*0.5*op,I*0.5*op]
    jednacina = np.linalg.solve(sistem,resenje)
    jednacina1 = jednacina[0]/(I*op)
#    jednacina2 = jednacina[1]/(I*op)
    return(jednacina1)

def integral(func,a,b,x):
    def realnaFja(x):
        return scipy.real(func(x))
    def imaginarnaFja(x):
        return scipy.imag(func(x))
    integralRe = integrate.quad(realnaFja,a,b)
    integralIm = integrate.quad(imaginarnaFja,a,b)
    return (integralRe[0]+I*integralIm[0],integralRe[1]+I*integralIm[1])

Ep = integral(f,-np.inf,np.inf,w)
print(Ep)


