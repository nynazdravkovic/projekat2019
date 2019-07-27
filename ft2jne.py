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
oc=1.
op=0.5
W = 0.01
gamma01=1.
gamma02=1.
I = complex(0,1)
n=10
c=1
udaljenost = np.linspace(0,5,n)
vreme = np.linspace(0,5,n)
E = np.empty((n,n))
epsilon = 1

#vraca ro02 sa tildom
def f (w):
    sistem = [[I*w+(-gamma01+I*(dc-dp)),-I*0.5*oc],[-I*0.5*oc,I*w+(-gamma02 - I*dp)]]
    resenje = [0,I*0.5*op]
    jednacina = np.linalg.solve(sistem,resenje)
#    jednacina1 = jednacina[0]/(I*op)
    jednacina2 = jednacina[1]/(I*op)
    return(jednacina2)

#funkcija koja racuna povrsinu ispod funkcije za odvojei imaginarni i realni deo
def integral(func,a,b):
    def realnaFja(x):
        return scipy.real(func(x))
    def imaginarnaFja(x):
        return scipy.imag(func(x))
    integralRe = integrate.quad(realnaFja,a,b)
    integralIm = integrate.quad(imaginarnaFja,a,b)
    return (integralRe[0]+I*integralIm[0],integralRe[1]+I*integralIm[1])

#ovo smo izracunale rucno 
def gaus(w,W):
    a = epsilon*W/np.sqrt(2*np.log(2))*np.exp(-w**2*W**2/(8*np.log(8)))/np.sqrt(2)
    return(a)

#Ep0 je niz koji se odmah racuna, to je gaus koji je prethodni prosao kroz fft 
#da bi bio gaus treba da t ide u neg vrednosti 
#k je mesto u nizu a t je vreme zato sto je vreme i negaivno a niz je izracunat tim vremenom pre
def ispodIntegrala(w,t,k, z):
    return(gaus(w,W)*np.exp(-I*w*(t-z/c)-f(w)*z))
    
def resavanje(distance,time,W):
    matricaResenjaRe = np.zeros((len(time),len(distance)))
    matricaResenjaIm = np.zeros((len(time),len(distance)))
    matricaModuo=np.zeros((len(time),len(distance)))
    for  i in range(len(time)):
        t = time[i]
        for j in range(len(distance)):
            z = distance[j]
            Ep = integral(lambda w: ispodIntegrala(w,t,i,z),-np.inf,np.inf)
            matricaResenjaRe[i][j] = 1/(np.sqrt(2*np.pi))*Ep[0].real
            matricaResenjaIm[i][j] = 1/(np.sqrt(2*np.pi))*Ep[0].imag
            matricaModuo[i][j]=np.sqrt(matricaResenjaRe[i][j]**2+matricaResenjaIm[i][j]**2)
    return(matricaResenjaRe,matricaResenjaIm,matricaModuo)
    

poljeRe = resavanje(udaljenost,vreme,W)[0]
poljeIm =  resavanje(udaljenost,vreme,W)[1]
poljemoduo=resavanje(udaljenost,vreme,W)[2]
#np.savetxt('Ep1.csv',polje, delimiter=',')
N=100
T = np.meshgrid(vreme,udaljenost)
#plt.contourf(T[1],T[0],poljeRe)
plt.contourf(T[1], T[0], poljeRe, N, cmap="viridis")
plt.colorbar()
plt.show()
plt.contourf(T[1], T[0], poljeIm, N, cmap="viridis")
plt.colorbar()
plt.show()
plt.contourf(T[1], T[0], poljemoduo, N, cmap="viridis")
plt.colorbar()
plt.show()
#

