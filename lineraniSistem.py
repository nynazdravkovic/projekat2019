# -*- coding: utf-8 -*-
"""
Created on Sat Jul 27 19:45:13 2019

@author: nina
"""
import numpy as np
import matplotlib.pyplot as plt

#master jednacine izvedene iz hamiltonijana za sistem od 3 nivoa i 3 lasera 
##resavanje sistema linearnih jednacina u redu 
I = complex(0,1)
Gamma01 = 1.
Gamma02 = 1.
Gamma12 = 1.
gamma01 = 1.
gamma02 = 1.
gamma10 = 1.
gamma12 = 1.
gamma20 = 1.
gamma21 = 1.
dc1 = 0.1
dc2 = 0.1
#dp = 0.
oc1 = 10.
oc2 = 5.
op = 0.5
#dp = 0.1

#f1 = Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*oc1*ro01 - 0.5*oc1*ro10)
#f2 = -gamma01*ro01 - 1.0*I*(-dc1*ro01 + 0.5*oc1*ro00 - 0.5*oc1*ro11 + 0.5*oc2*ro02)
#f3 = -gamma02*ro02 - 1.0*I*(dp*ro02 - 0.5*oc1*ro12 + 0.5*oc2*ro01)
#f4 = -gamma10*ro10 - 1.0*I*(dc1*ro10 - 0.5*oc1*ro00 + 0.5*oc1*ro11 - 0.5*oc2*ro20)
#f5 = -Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(-0.5*oc1*ro01 + 0.5*oc1*ro10 + 0.5*oc2*ro12 - 0.5*oc2*ro21)
#f6 = -gamma12*ro12 - 1.0*I*(dc1*ro12 + dp*ro12 - 0.5*oc1*ro02 + 0.5*oc2*ro11 - 0.5*oc2*ro22)
#f7 = -gamma20*ro20 - 1.0*I*(-dp*ro20 + 0.5*oc1*ro21 - 0.5*oc2*ro10)
#f8 = -gamma21*ro21 - 1.0*I*(-dc1*ro21 - dp*ro21 + 0.5*oc1*ro20 - 0.5*oc2*ro11 + 0.5*oc2*ro22)
#f9 = -Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc2*ro12 + 0.5*oc2*ro21)
rhoIm=[]
rhoRe=[]
rho00=[]
rho11=[]
rho22=[]
delta = np.linspace(-20,20,400)
zbirDijagonale = []
for i in range(400):
    dp = delta[i]
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
# ovde resavam jednacine za prvi stepen a ono sto stoji uz op uzimam iz 
#prethodne matrice jer uz njega stoji nesto sto mnoz lambda na 0 :)
    prviStepen=[[1,0,0,0,1,0,0,0,1],
                 [-I*oc1/2,-gamma01+dc1,-I*oc2/2,0,I*oc1/2,0,0,0,0],
                 [0,-I*oc2/2,-gamma02-I*dp,0,0,I*oc1/2,0,0,0],
                 [I*oc1/2,0,0,-gamma10-I*dc1,-I*oc1/2,0,I*oc2/2,0,0],
                 [0,I*oc1/2,0,-I*oc1/2,-Gamma01,-I*oc2/2,0,I*oc2/2,Gamma12],
                 [0,0,I*oc1/2,0,-I*oc2/2,-gamma12-I*(dc1+dp),0,0,I*oc2/2],
                 [0,0,0,I*oc2/2,0,0,-gamma20+I*dp,-I*oc1/2,0],
                 [0,0,0,0,I*oc2/2,0,-I*oc1/2,-gamma21+I*(dc1+dp),-I*oc2/2],
                 [0,0,0,0,0,I*oc2/2,0,-I*oc2/2,-Gamma02-Gamma12]]
    resenjePrvog = [1,
                    I*op*0.5*resenje1[7],
                    I*op*0.5*(resenje1[0]-resenje1[8]),
                    I*op*0.5*resenje1[5],
                    0,
                    I*op*0.5*resenje1[3],
                    I*op*0.5*(resenje1[8]-resenje1[0]),
                    I*op*0.5*resenje1[1],
                    I*op*0.5*(resenje1[6]-resenje1[2])]
    resenje2 = np.linalg.solve(prviStepen,resenjePrvog)
    zbirDijagonale.append(resenje2[0]+resenje2[4]+resenje2[8])
    rhoRe.append(resenje2[6].real)
    rhoIm.append(resenje2[6].imag)
    rho00.append(resenje2[0].real)
    rho11.append(resenje2[4].real)
    rho22.append(resenje2[8].real)
    
plt.plot(delta,rhoRe)
plt.plot(delta,rhoIm)
plt.show()
plt.plot(delta,rho00)
plt.plot(delta,rho11)
plt.plot(delta,rho22)
plt.plot(delta,zbirDijagonale)
plt.show()
