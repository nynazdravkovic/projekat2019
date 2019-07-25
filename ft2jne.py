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

dc=0.1
dp=0.1
oc=0.5
op=0.5
W = 1
gamma01=1.
gamma02=1.

I = complex(0,1)
n=100
c=1
udaljenost = np.linspace(0,10,n,dtype = 'int16')
vreme = np.linspace(-5,5,n,dtype = 'int16')
E = np.empty((n,n))
epsilon = 1

#vraca ro02 sa tildom
def f (w):
    sistem = [[I*w+(-gamma01+I*(dc-dp)),-I*0.5*oc],[-I*0.5*oc,I*w+(-gamma02 - I*dp)]]
    resenje = [-I*0.5*op,I*0.5*op]
    jednacina = np.linalg.solve(sistem,resenje)
#    jednacina1 = jednacina[0]/(I*op)
    jednacina2 = jednacina[1]/(I*op)
    return(np.conj(jednacina2))

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
#da bi bio gaus treba da t ide u neg vrednosti 
#k je mesto u nizu a t je vreme zato sto je vreme i negaivno a niz je izracunat tim vremenom pre
def ispodIntegrala(w, Ep0,t,k, z):
    return(Ep0[k]*np.exp(-I*w*(t-z/c)-f(w)*z))
    
def resavanje(distance,time,W):
    matricaResenja = np.zeros((len(time),len(distance)),dtype = 'complex')
    Ep0 = gaus(len(time),W)
    for  i in range(len(time)):
        t = time[i]
        for j in range(len(distance)):
            z = distance[j]
            Ep = integral(lambda w: ispodIntegrala(w,Ep0,t,i,z),-np.inf,np.inf)
            matricaResenja[i][j] = 1/(2*np.pi)*Ep[0]
    return(matricaResenja)
   
polje = resavanje(vreme,udaljenost,W)
#np.savetxt('Ep1.csv',polje, delimiter=',')

T = np.meshgrid(vreme,udaljenost)
plt.contourf(T[0],T[1],polje)
plt.colorbar()
plt.show()
#

