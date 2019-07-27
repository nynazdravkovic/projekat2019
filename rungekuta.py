# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 11:27:34 2019
@author: nina
"""

#ro'00 = ro11 - I*(o*ro01/2 - o*ro10/2)
#ro'01 = -ro01 - I*(-d*ro01 + o*ro00/2 - o*ro11/2)
#ro'10 = -ro10 - I*(d*ro10 - o*ro00/2 + o*ro11/2)
#ro'11 = -ro11 - I*(-o*ro01/2 + o*ro10/2)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
Gamma02=1.
Gamma12=1.
Gamma01=1.
gamma01=1.
gamma02=1.
gamma10=1.
gamma12=1.
gamma20=1.
gamma21=1.

ro00 = 1.
ro01 = 0.
ro02 = 0.
ro10 = 0.
ro11 = 0.
ro12 = 0.
ro20 = 0.
ro21 = 0.
ro22 = 0.

op = 0.5
oc = 0.5
dp = 0.
dc = 0.
oc1 = 0.1
oc2 = 0.1
dc1 = 0.
dc2 = 0.

I = complex(0,1)
pocetniUslovi = [1.,0.,0.,0.,0.,0.,0.,0.,0.]
matricaGustineIm = []
matricaGustineRe = []
t = 0.
h = 0.02
ro00re = []
time = []
#ako se koristi ova funkcija koriste se master jednacine za sistem sa 3 nivoa 
# i 3lasera
#def mnozenjeK(ro):
#    k1 = Gamma01*ro[4]+ Gamma02*ro[8] - 1.0*I*(0.5*op*ro[2] - 0.5*op*ro[6])
#    k2 = -gamma01*ro[1] - 1.0*I*(0.5*oc*ro[2] - 0.5*op*ro[7] - ro[1]*(dc - dp))
#    k3 = -gamma02*ro[2] - 1.0*I*(dp*ro[2] + 0.5*oc*ro[1] + 0.5*op*ro[0] - 0.5*op*ro[8])
#    k4 = -gamma10*ro[3] - 1.0*I*(-0.5*oc*ro[6] + 0.5*op*ro[5] + ro[3]*(dc - dp)) 
#    k5 = -Gamma01*ro[4] + Gamma12*ro[8] - 1.0*I*(0.5*oc*ro[5] - 0.5*oc*ro[7])
#    k6= -gamma12*ro[5] - 1.0*I*(dp*ro[5] + 0.5*oc*ro[4] - 0.5*oc*ro[8] + 0.5*op*ro[3] + ro[5]*(dc - dp))
#    k7 = -gamma20*ro[6] - 1.0*I*(-dp*ro[6] - 0.5*oc*ro[3] - 0.5*op*ro[0] + 0.5*op*ro[8])
#    k8 = -gamma21*ro[7] - 1.0*I*(-dp*ro[7] - 0.5*oc*ro[4] + 0.5*oc*ro[8] - 0.5*op*ro[5] - ro[7]*(dc - dp))
#    k9 = -Gamma02*ro[8] - Gamma12*ro[8] - 1.0*I*(-0.5*oc*ro[5] + 0.5*oc*ro[7] - 0.5*op*ro[2] + 0.5*op*ro[6])
#    return(k1,k2,k3,k4,k5,k6,k7,k8,k9)
#ako se koristi ova funkcija koriste se master jednacine za sistem sa 3 nivao i 
# 3 lasera 
def mnozenjeK(ro):        
    f1 = Gamma01*ro[4] + Gamma02*ro[8] - 1.0*I*(0.5*oc1*ro[1] - 0.5*oc1*ro[3])
    f2 = -gamma01*ro[1] - 1.0*I*(-dc1*ro[1] + 0.5*oc1*ro[0] - 0.5*oc1*ro[4] + 0.5*oc2*ro[2])
    f3 = -gamma02*ro[2] - 1.0*I*(dp*ro[2] - 0.5*oc1*ro[5] + 0.5*oc2*ro[1])
    f4 = -gamma10*ro[3] - 1.0*I*(dc1*ro[3] - 0.5*oc1*ro[0] + 0.5*oc1*ro[4] - 0.5*oc2*ro[6])
    f5 = -Gamma01*ro[4] + Gamma12*ro[8] - 1.0*I*(-0.5*oc1*ro[1] + 0.5*oc1*ro[3] + 0.5*oc2*ro[5] - 0.5*oc2*ro[7])
    f6 = -gamma12*ro[5] - 1.0*I*(dc1*ro[6] + dp*ro[5] - 0.5*oc1*ro[2] + 0.5*oc2*ro[4] - 0.5*oc2*ro[8])
    f7 = -gamma20*ro[6] - 1.0*I*(-dp*ro[6] + 0.5*oc1*ro[7] - 0.5*oc2*ro[3])
    f8 = -gamma21*ro[7] - 1.0*I*(-dc1*ro[7] - dp*ro[7] + 0.5*oc1*ro[6] - 0.5*oc2*ro[4] + 0.5*oc2*ro[8])
    f9 = -Gamma02*ro[2] - Gamma12*ro[8] - 1.0*I*(-0.5*oc2*ro[5] + 0.5*oc2*ro[7])
    return(f1,f2,f3,f4,f5,f6,f7,f8,f9)
def dodavanjeRo(k,ro):
    ro00 = ro[0]+(h/2.)+k[0]
    ro01 = ro[1]+(h/2.)+k[1]
    ro02 = ro[2]+(h/2.)+k[2]
    ro10 = ro[3]+(h/2.)+k[3]
    ro11 = ro[4]+(h/2.)+k[4]
    ro12 = ro[5]+(h/2.)+k[5]
    ro20 = ro[6]+(h/2.)+k[6]
    ro21 = ro[7]+(h/2.)+k[7]
    ro22 = ro[8]+(h/2.)+k[8]
    return(ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22)
rho11 = []
rho22 = []
rho00 = []
rho02re = []
rho02im = []
while (t<4):
    K1 = mnozenjeK(pocetniUslovi)
    RO1 = dodavanjeRo(K1,pocetniUslovi)
    K2 = mnozenjeK(RO1)
    RO2 = dodavanjeRo(K2,pocetniUslovi)
    K3 = mnozenjeK(RO2)
    RO3 = dodavanjeRo(K3,pocetniUslovi)
    K4 = mnozenjeK(RO3)
    t = t + (h/2.)  
    for j in range(9):
        pocetniUslovi[j] += (h/6.)*(K1[j]+(2.*K2[j])+(2.*K3[j])+K4[j])
    time.append(t)  
    rho00.append(pocetniUslovi[0].real)
    rho11.append(pocetniUslovi[4].real)
    rho22.append(pocetniUslovi[8].real)
    rho02re.append(pocetniUslovi[2].real)
    rho02im.append(pocetniUslovi[2].imag)
plt.plot(time,rho00, label = 'rho00')
plt.plot(time,rho11, label = 'rho11')
plt.plot(time,rho22, label = 'rho22')
plt.legend()
plt.xlabel("vreme")
plt.ylabel("rho")
plt.show()

plt.plot(time,rho02re, label = 'realni deo')
plt.plot(time,rho02im, label = 'imaginarni deo')
plt.legend()
plt.xlabel("vreme")
plt.ylabel("rho02")
plt.show()