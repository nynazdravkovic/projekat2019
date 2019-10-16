# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 21:38:32 2019

@author: nina
"""


import numpy as np
import matplotlib.pyplot as plt
import scipy
import scipy.integrate as integrate
from scipy.fftpack import fft, ifft
#W konfiguracija
dc2=1.
dc1=1.
dp=1.
oc2=10.
oc1=10.
op=0.5
W = 4.
Gamma01 = 1.
Gamma02 = 1.
Gamma12 = 1.
Gamma03 = 1.
Gamma13 = 1.
Gamma23 = 1.
gamma23 = 1.
gamma01 = 1.
gamma02 = 1.
gamma10 = 1.
gamma12 = 1.
gamma20 = 1.
gamma21 = 1.
gamma03 = 1.
gamma13 = 1.
gamma13 = 1.
gamma30 = 1.
gamma31 = 1.
gamma32 = 1.
gamma33 = 1.

I = complex(0,1)
n=10
c=1
udaljenost = np.linspace(0,5,n)
vreme = np.linspace(0,5,n)
E = np.empty((n,n))
epsilon = 1

#vraca ro02 sa tildom
def f (w):
    nultiStepen=[[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
                 [-I*oc1,-gamma01+I*dc1+I*w,0,0,0,I*oc1,0,0,I*oc2,0,0,0,0,0,0,0],
                 [-I*oc2,0,-gamma02+I*dc2+I*w,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0,0],
                 [0,0,0,I*dp-gamma03+I*w,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0],
                 [I*oc1,0,0,0,-gamma10-I*dc1+I*w,-I*oc1,-oc2*I,0,0,0,0,0,0,0,0,0],
                 [0,I*oc1,0,0,-I*oc1,-Gamma01+I*w,0,0,0,0,Gamma12,Gamma13,0,0,0,0],
                 [0,0,I*oc1,0,-I*oc2,0,dc2-dc1-gamma12+I*w,0,0,0,0,0,0,0,0,0],
                 [0,0,0,I*oc1,0,0,0,dp-dc1-gamma13+I*w,0,0,0,0,0,0,0,0],
                 [-I*oc2,0,0,0,0,0,0,0,-gamma20-I*dc2+I*w,-I*oc1,-I*oc2,0,0,0,0,0],
                 [0,I*oc2,0,0,0,0,0,0,-I*oc1,-gamma21+I*(dc1-dc2)+I*w,0,0,0,0,0,0],
                 [0,0,I*oc2,0,0,0,0,0,-I*oc2,0,-Gamma02-Gamma12+I*w,0,0,0,0,Gamma23],
                 [0,0,0,oc2*I,0,0,0,0,0,0,0,-gamma23+I*(dp-dc2),0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-gamma30+I*w-I*dp,-I*oc1,-I*oc2,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc1,-gamma31-I*(dp-dc1)+I*w,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc2,0,-gamma32-I*(dp-dc2)+I*w,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Gamma03]]
                    
    
    nule = np.zeros(16,dtype='complex')
    nule[0] = 1
    resenje1 = np.linalg.solve(nultiStepen,nule)
# ovde resavam jednacine za prvi stepen a ono sto stoji uz op uzimam iz 
#prethodne matrice jer uz njega stoji nesto sto mnoz lambda na 0 :)
    prviStepen=[[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
                 [-I*oc1,-gamma01+I*dc1+I*w,0,0,0,I*oc1,0,0,I*oc2,0,0,0,0,0,0,0],
                 [-I*oc2,0,-gamma02+I*dc2+I*w,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0,0],
                 [0,0,0,I*dp-gamma03+I*w,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0],
                 [I*oc1,0,0,0,-gamma10-I*dc1+I*w,-I*oc1,-oc2*I,0,0,0,0,0,0,0,0,0],
                 [0,I*oc1,0,0,-I*oc1,-Gamma01+I*w,0,0,0,0,Gamma12,Gamma13,0,0,0,0],
                 [0,0,I*oc1,0,-I*oc2,0,dc2-dc1-gamma12+I*w,0,0,0,0,0,0,0,0,0],
                 [0,0,0,I*oc1,0,0,0,dp-dc1-gamma13+I*w,0,0,0,0,0,0,0,0],
                 [-I*oc2,0,0,0,0,0,0,0,-gamma20-I*dc2+I*w,-I*oc1,-I*oc2,0,0,0,0,0],
                 [0,I*oc2,0,0,0,0,0,0,-I*oc1,-gamma21+I*(dc1-dc2)+I*w,0,0,0,0,0,0],
                 [0,0,I*oc2,0,0,0,0,0,-I*oc2,0,-Gamma02-Gamma12+I*w,0,0,0,0,Gamma23],
                 [0,0,0,oc2*I,0,0,0,0,0,0,0,-gamma23+I*(dp-dc2),0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-gamma30+I*w-I*dp,-I*oc1,-I*oc2,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc1,-gamma31-I*(dp-dc1)+I*w,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc2,0,-gamma32-I*(dp-dc2)+I*w,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Gamma03]]
    resenjePrvog = [1,
                    -I*op*resenje1[12],
                    -I*op*resenje1[13],
                    I*op*(resenje1[0]-resenje1[10]),
                    -I*op*resenje1[7],
                    0,
                    0,
                    -I*op*resenje1[4],
                    -I*op*resenje1[11],
                    0,
                    0,
                    I*op*resenje1[8],
                    I*op*(resenje1[15]-resenje1[0]),
                    I*op*resenje1[1],
                    I*op*resenje1[2],
                    I*op*(resenje1[4]-resenje1[12])]
    resenje2 = np.linalg.solve(prviStepen,resenjePrvog)
    return(resenje2[12]/(I*op))

#funkcija koja racuna povrsinu ispod funkcije za odvojei imaginarni i realni deo
def integral(func,a,b):
    def realnaFja(x):
        return scipy.real(func(x))
    def imaginarnaFja(x):
        return scipy.imag(func(x))
    integralRe = integrate.quad(realnaFja,a,b)
    integralIm = integrate.quad(imaginarnaFja,a,b)
    return (integralRe[0]+I*integralIm[0])

#ovo smo izracunale rucno 
def gaus(w,W):
    a = epsilon*W/np.sqrt(2*np.log(2))*np.exp(-(w)**2*(W)**2/(8*np.log(8)))/np.sqrt(2)
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
#            print(Ep)
            matricaResenjaRe[i][j] = 1/(np.sqrt(2*np.pi))*Ep.real
            matricaResenjaIm[i][j] = 1/(np.sqrt(2*np.pi))*Ep.imag
            matricaModuo[i][j]=np.sqrt(matricaResenjaRe[i][j]**2+matricaResenjaIm[i][j]**2)
    return(matricaResenjaRe,matricaResenjaIm,matricaModuo)
    
polje = resavanje(udaljenost,vreme,W)
poljeIm = polje[1]
poljeRe = polje[0]
poljeModuo = polje[2]

#np.savetxt('Ep1.csv',polje, delimiter=',')
N=n*n
T = np.meshgrid(vreme,udaljenost)
#plt.contourf(T[1],T[0],poljeRe)
plt.contourf(T[1], T[0], poljeRe, N, cmap="viridis")
plt.colorbar()

plt.show()
plt.contourf(T[1], T[0], poljeIm, N, cmap="viridis")
plt.colorbar()
plt.show()
plt.contourf(T[1], T[0], poljeModuo, N, cmap="viridis")
plt.colorbar()
plt.title("Moduo probnog polja, Ep")
plt.xlabel("x[m]")
plt.ylabel("t[s]")
plt.show()
#
#
