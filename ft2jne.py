# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 11:01:39 2019

@author: nina
"""

import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.integrate as integrate
from scipy.fftpack import fft, ifft

#za resavanje simbolicki 
#import sympy
#from sympy import *
#from sympy.solvers.solvers import solve_linear_system
#dc=sympy.Symbol('dc')
#dp=sympy.Symbol('dp')
#oc=sympy.Symbol('oc')
#op=sympy.Symbol('op')
#ro00=sympy.Symbol('ro00')
#ro01=sympy.Symbol('ro01')
#ro02=sympy.Symbol('ro02')
#gamma01=sympy.Symbol('gamma01')
#gamma02=sympy.Symbol('gamma02')
#I = sympy.I
#za direktno resavanje
dc=0.
dp=0.
oc=0.5
op=0.5
W = 0.1
gamma01=1.
gamma02=1.
#w = 1.
I = complex(0,1)
n=100
omega = np.linspace(0,100,100)
c=1
udaljenost = []
vreme = []
E = np.empty((100,100))
epsilon = 1
#jednacina = sympy.solve([I*w*ro01-I*0.5*oc*ro02+I*0.5*op+ro01*(-gamma01+I*(dc-dp)),I*w*ro02-I*0.5*oc*ro01-I*0.5*op+ro02*(-gamma02 - I*dp)],[ro01,ro02])
#print(jednacina[ro01])
def f (w):
    sistem = [[I*w+(-gamma01+I*(dc-dp)),-I*0.5*oc],[-I*0.5*oc,I*w+(-gamma02 - I*dp)]]
    resenje = [-I*0.5*op,I*0.5*op]
    jednacina = np.linalg.solve(sistem,resenje)
#    jednacina1 = jednacina[0]/(I*op)
    jednacina2 = jednacina[1]/(I*op)
    return(jednacina2)
print()
#funkcija koja racuna povrsinu ispod funkcije za odvojei imaginarni i realni deo
def integral(func,a,b):
    def realnaFja(x):
        return scipy.real(func(x))
    def imaginarnaFja(x):
        return scipy.imag(func(x))
    integralRe = integrate.quad(realnaFja,a,b)
    integralIm = integrate.quad(imaginarnaFja,a,b)
    return (integralRe[0]+I*integralIm[0],integralRe[1]+I*integralIm[1])

#vraca mi niz E(0,t) koji je prosao kroz fft sigurno radi
def gaus(time,W):
    niz = []
    for i in range (time):
        t = i
        a = epsilon*np.exp(-2*np.log(2)*(t**2/W**2))    
        niz.append(a)
    a = fft(niz)
    return(a)
#Ep0 je niz koji se odmah racuna, to je gaus koji je prethodni prosao kroz fft 
def ispodIntegrala(w, Ep0,t, z):
    return(Ep0[t]*np.exp(-I*w*(t-z/c)-np.conj(f(w))*z))
def resavanje(distance,time,W):
    matricaResenja = np.zeros((time,distance),dtype = 'complex')
    for  i in range(time):
        t = i
        Ep0 = gaus(time,W)
        for j in range(distance):
            z = j
            Ep = integral(lambda w: ispodIntegrala(w,Ep0,t,z),0,np.inf)
            matricaResenja[i][j] = Ep[0]
    return(matricaResenja)
   

#def gaus(time,W):
#    a = epsilon*np.exp(-2*np.log(time**2/W**2))    
#    return(a)
##Ep0 je niz koji se odmah racuna, to je gaus koji je prethodni prosao kroz fft 
#def ispodIntegrala(w,E0,t,z):
#    return(E0*np.exp(-I*w*(t-z/c)-np.conj(f(w))*z))
#def resavanje(distance,time,W):
#   matricaResenja = np.zeros((time,distance),dtype = 'complex')
#   matricaGreski = np.zeros((time,distance),dtype = 'complex')
#   for  i in range(time):
#       t = i
#       vreme.append(t)
#       udaljenost.append(t)
#       Ep0 = gaus(time,W)
#       for j in range(distance):
#           z = j
#           Ep = integral(lambda w: ispodIntegrala(w,Ep0,t,z),-np.inf,np.inf)
#           matricaResenja[i][j] = Ep[0]
#           matricaGreski[i][j] = Ep[1]
#   return(matricaResenja,matricaGreski)

polje = resavanje(n,n,W)
poljeGreske = resavanje(n,n,W)
#np.savetxt('Ep1.csv',polje, delimiter=',')
T,D = np.meshgrid(vreme,udaljenost)
plt.contourf(polje)
plt.show()


