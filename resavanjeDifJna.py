# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:23:48 2019

@author: nina
"""
import numpy as np
from scipy.integrate import odeint, solve_ivp
import matplotlib.pyplot as plt


I=complex(0,1)
op=0.5
oc=0.5
dp=0.
dc=0.

gamma00=1.
gamma01=1.
gamma02=1.
gamma10=1.
gamma11=1.
gamma12=1.
gamma20=1.
gamma21=1.
gamma22=1.

Gamma00=1.
Gamma01=1.
Gamma02=1.
Gamma10=1.
Gamma11=1.
Gamma12=1.
Gamma20=1.
Gamma21=1.
Gamma22=1.
I=complex(0,1)
op=0.5
oc1=0.5
oc2 = 0.5
dp=0.
dc1=0.
dc2 = 0.
#    prvo resavam jednacine za nulti stepen ovde
nultiStepen=[[1,0,0,0,1,0,0,0,1],
             [-I*oc1/2,-gamma01+dc1,-I*oc2/2,0,I*oc1/2,0,0,0,0],
             [0,-I*oc2/2,-gamma02-I*dp,0,0,I*oc1/2,0,0,0],
             [I*oc1/2,0,0,-gamma10-I*dc1,-I*oc1/2,0,I*oc2/2,0,0],
             [0,I*oc1/2,0,-I*oc1/2,-Gamma01,-I*oc2/2,0,I*oc2/2,Gamma12],
             [0,0,I*oc1/2,0,-I*oc2/2,-gamma12-I*(dc1+dp),0,0,I*oc2/2],
             [0,0,0,I*oc2/2,0,0,-gamma20+I*dp,-I*oc1/2,0],
             [0,0,0,0,I*oc2/2,0,-I*oc1/2,-gamma21+I*(dc1+dp),-I*oc2/2],
             [0,0,0,0,0,I*oc2/2,0,-I*oc2/2,-Gamma02-Gamma12]]

nule = np.zeros(9)
nule[0] = 1
resenje1 = np.linalg.solve(nultiStepen,nule)

def model(t,z):
    ro00 = z[0]
    ro01 = z[1]
    ro02 = z[2]
    ro10 = z[3]  
    ro11 = z[4]
    ro12 = z[5]
    ro20 = z[6]
    ro21 = z[7] 
    ro22 = z[8]   
    dro00dt = Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*oc1*ro01 - 0.5*oc1*ro10 + 0.5*op*ro02 - 0.5*op*ro20)
    dro01dt = -gamma01*ro01 - 1.0*I*(-dc1*ro01 + 0.5*oc1*ro00 - 0.5*oc1*ro11 + 0.5*oc2*ro02 - 0.5*op*ro21)
    dro02dt = -gamma02*ro02 - 1.0*I*(dp*ro02 - 0.5*oc1*ro12 + 0.5*oc2*ro01 + 0.5*op*ro00 - 0.5*op*ro22)
    dro10dt = -gamma10*ro10 - 1.0*I*(dc1*ro10 - 0.5*oc1*ro00 + 0.5*oc1*ro11 - 0.5*oc2*ro20 + 0.5*op*ro12)
    dro11dt = -Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(-0.5*oc1*ro01 + 0.5*oc1*ro10 + 0.5*oc2*ro12 - 0.5*oc2*ro21)
    dro12dt = -gamma12*ro12 - 1.0*I*(dc1*ro12 + dp*ro12 - 0.5*oc1*ro02 + 0.5*oc2*ro11 - 0.5*oc2*ro22 + 0.5*op*ro10)
    dro20dt = -gamma20*ro20 - 1.0*I*(-dp*ro20 + 0.5*oc1*ro21 - 0.5*oc2*ro10 - 0.5*op*ro00 + 0.5*op*ro22)
    dro21dt = -gamma21*ro21 - 1.0*I*(-dc1*ro21 - dp*ro21 + 0.5*oc1*ro20 - 0.5*oc2*ro11 + 0.5*oc2*ro22 - 0.5*op*ro01)
    dro22dt = -Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc2*ro12 + 0.5*oc2*ro21 - 0.5*op*ro02 + 0.5*op*ro20)
    dzdt = [dro00dt, dro01dt, dro02dt, dro10dt, dro11dt, dro12dt, dro20dt, dro21dt, dro22dt]
    return dzdt

z0 = np.array([1., 0., 0., 0., 0., 0., 0., 0., 0.], dtype='complex')
t = np.linspace(0,5,50000)
z = solve_ivp(model, (0., 5.), resenje1)
t=z.t
z=np.transpose(z.y)

ro00 = z[:,0]
ro01 = z[:,1]
ro02 = z[:,2]
ro10 = z[:,3]
ro11 = z[:,4]
ro12 = z[:,5]
ro20 = z[:,6]
ro21 = z[:,7]
ro22 = z[:,8]

plt.plot(t, ro00.real+ro11.real+ro22.real, label='zbir realnih delova ro22 i rho11 i rho00')
plt.ylabel('\rho00.real+\rho11.real')
plt.xlabel('vreme')
plt.legend(loc='best')
plt.show()
plt.plot(t,ro00.real,'g:',label='ro00 realni deo')
plt.plot(t,ro01.real,'b:',label='ro01 realni deo')
plt.plot(t,ro02.real,'r:',label='ro02 realni deo')
plt.plot(t,ro10.real,'m:',label='ro10 realni deo')
plt.plot(t,ro11.real,'g-',label='ro11 realni deo')
plt.plot(t,ro12.real,'b-',label='ro12 realni deo')
plt.plot(t,ro20.real,'r-',label='ro20 realni deo')
plt.plot(t,ro21.real,'m-',label='ro21 realni deo')
plt.plot(t,ro22.real,'m--',label='ro22 realni deo')
#plt.plot(t, ro00.real+ro11.real)
plt.ylabel('realni delovi \rhomn')
plt.xlabel('vreme')
plt.legend(loc='best')
plt.show()

plt.plot(t,ro00.imag,'g:',label='ro00 imaginarni deo')
plt.plot(t,ro01.imag,'b:',label='ro01 imaginarni deo')
plt.plot(t,ro02.imag,'r:',label='ro02 imaginarni deo')
plt.plot(t,ro10.imag,'m:',label='ro10 imaginarni deo')
plt.plot(t,ro11.imag,'g-',label='ro11 imaginarni deo')
plt.plot(t,ro12.imag,'b-',label='ro12 imaginarni deo')
plt.plot(t,ro20.imag,'r-',label='ro20 imaginarni deo')
plt.plot(t,ro21.imag,'m-',label='ro21 imaginarni deo')
plt.plot(t,ro22.imag,'m--',label='ro22 imaginarni deo')
#plt.plot(t, ro00.real+ro11.real)
plt.ylabel('imaginarni delovi rhomn')
plt.xlabel('vreme')
plt.legend(loc='best')
plt.show()

#jedan ro00 i ro11, posle im i real
