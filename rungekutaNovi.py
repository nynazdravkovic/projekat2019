# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 23:27:54 2019

@author: nina
"""


#ro'00 = ro11 - I*(o*ro01/2 - o*ro10/2)
#ro'01 = -ro01 - I*(-d*ro01 + o*ro00/2 - o*ro11/2)
#ro'10 = -ro10 - I*(d*ro10 - o*ro00/2 + o*ro11/2)
#ro'11 = -ro11 - I*(-o*ro01/2 + o*ro10/2)
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
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
op = 0.5
oc = 0.5
dp = 0.
dc = 0.
oc1 = 0.5
oc2 = 0.5
dc1 = 0.
dc2 = 0.

I = complex(0,1)

matricaGustineIm = []
matricaGustineRe = []
rho00re = []
rho01re = []
rho02re = []
rho10re = []
rho11re = []
rho12re = []
rho20re = []
rho21re = []
rho22re = []

rho00im = []
rho01im = []
rho02im = []
rho10im = []
rho11im = []
rho12im = []
rho20im = []
rho21im = []
rho22im = []
t = 0.
h = 0.02
ro00re = []
time = []
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
pocetniUslovi = np.linalg.solve(nultiStepen,nule)
ro00 = pocetniUslovi[0]
ro01 = pocetniUslovi[1]
ro02 = pocetniUslovi[2]
ro10 = pocetniUslovi[3]
ro11 = pocetniUslovi[4]
ro12 = pocetniUslovi[5]
ro20 = pocetniUslovi[6]
ro21 = pocetniUslovi[7]
ro22 = pocetniUslovi[8]
#ako se koristi ova funkcija koriste se master jednacine za sistem sa 3 nivoa 
# i 3lasera
#pocetniUslovi = np.zeros(9,dtype = 'complex')
#pocetniUslovi[0] = 1
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
## 3 lasera 
def mnozenjeK(ro):        
    f1 = Gamma01*ro[4] + Gamma02*ro[8] - 1.0*I*(0.5*oc1*ro[1] - 0.5*oc1*ro[3] + 0.5*op*ro[2] - 0.5*op*ro[6])
    f2 = -gamma01*ro[1] - 1.0*I*(-dc1*ro[1] + 0.5*oc1*ro[0] - 0.5*oc1*ro[4] + 0.5*oc2*ro[2] - 0.5*op*ro[7])
    f3 = -gamma02*ro[2] - 1.0*I*(dp*ro[2] - 0.5*oc1*ro[5] + 0.5*oc2*ro[1] + 0.5*op*ro[0] - 0.5*op*ro[8])
    f4 = -gamma10*ro[3] - 1.0*I*(dc1*ro[3] - 0.5*oc1*ro[0] + 0.5*oc1*ro[4] - 0.5*oc2*ro[6] + 0.5*op*ro[5])
    f5 = -Gamma01*ro[4] + Gamma12*ro[8] - 1.0*I*(-0.5*oc1*ro[1] + 0.5*oc1*ro[3] + 0.5*oc2*ro[5] - 0.5*oc2*ro[7])
    f6 = -gamma12*ro[5] - 1.0*I*(dc1*ro[5] + dp*ro[5] - 0.5*oc1*ro[2] + 0.5*oc2*ro[4] - 0.5*oc2*ro[8] + 0.5*op*ro[3])
  
    f7 = -gamma20*ro[6] - 1.0*I*(-dp*ro[6] + 0.5*oc1*ro[7] - 0.5*oc2*ro[3] - 0.5*op*ro[0] + 0.5*op*ro[8])
    f8 = -gamma21*ro[7] - 1.0*I*(-dc1*ro[7] - dp*ro[7] + 0.5*oc1*ro[6] - 0.5*oc2*ro[4] + 0.5*oc2*ro[8] - 0.5*op*ro[1])
    f9 = -Gamma02*ro[8] - Gamma12*ro[8] - 1.0*I*(-0.5*oc2*ro[5] + 0.5*oc2*ro[7] - 0.5*op*ro[2] + 0.5*op*ro[6])
    return([f1,f2,f3,f4,f5,f6,f7,f8,f9])

while (t<4):
    K1=np.array(mnozenjeK(pocetniUslovi))*h
    K2=np.array(mnozenjeK(pocetniUslovi+K1/2))*h
    K3=np.array(mnozenjeK(pocetniUslovi+K2/2))*h
    K4=np.array(mnozenjeK(pocetniUslovi+K3))*h
    t = t + (h)  
    for j in range(9):
        pocetniUslovi[j] += (1/6.)*(K1[j]+(2.*K2[j])+(2.*K3[j])+K4[j])
    time.append(t)  
    rho00re.append(pocetniUslovi[0].real)
    rho01re.append(pocetniUslovi[1].real)
    rho02re.append(pocetniUslovi[2].real)
    rho10re.append(pocetniUslovi[3].real)
    rho11re.append(pocetniUslovi[4].real)
    rho12re.append(pocetniUslovi[5].real)
    rho20re.append(pocetniUslovi[6].real)
    rho21re.append(pocetniUslovi[7].real)
    rho22re.append(pocetniUslovi[8].real)
    
    rho00im.append(pocetniUslovi[0].imag)
    rho01im.append(pocetniUslovi[1].imag)
    rho02im.append(pocetniUslovi[2].imag)
    rho10im.append(pocetniUslovi[3].imag)
    rho11im.append(pocetniUslovi[4].imag)
    rho12im.append(pocetniUslovi[5].imag)
    rho20im.append(pocetniUslovi[6].imag)
    rho21im.append(pocetniUslovi[7].imag)
    rho22im.append(pocetniUslovi[8].imag)

plt.plot(time,rho00re, 'g:', label = 'rho00')
plt.plot(time,rho01re,'b:' , label = 'rho01')
plt.plot(time,rho02re, 'r:', label = 'rho02')
plt.plot(time,rho10re,'m:' , label = 'rho10')
plt.plot(time,rho11re, 'g-', label = 'rho11')
plt.plot(time,rho12re, 'b-', label = 'rho12')
plt.plot(time,rho20re, 'r-', label = 'rho20')
plt.plot(time,rho21re, 'm', label = 'rho21')
plt.plot(time,rho22re, 'm--', label = 'rho22')
plt.legend()
#plt.title()
plt.xlabel("vreme")
plt.ylabel("rho")
plt.show()
plt.plot(time,rho00im, 'g:', label = 'rho00')
plt.plot(time,rho01im,'b:' , label = 'rho01')
plt.plot(time,rho02im, 'r:', label = 'rho02')
plt.plot(time,rho10im,'m:' , label = 'rho10')
plt.plot(time,rho11im, 'g-', label = 'rho11')
plt.plot(time,rho12im, 'b-', label = 'rho12')
plt.plot(time,rho20im, 'r-', label = 'rho20')
plt.plot(time,rho21im, 'm', label = 'rho21')
plt.plot(time,rho22im, 'm--', label = 'rho22')

plt.legend()
#plt.title("Grafik zavisnosti ")
plt.xlabel("vreme")
plt.ylabel("rho")
plt.show()
#
#plt.plot(time,rho02re, label = 'realni deo')
#plt.plot(time,rho02im, label = 'imaginarni deo')
#plt.legend()
#plt.xlabel("vreme")
#plt.ylabel("rho02")
#plt.show()
#plt.title("Grafik zavisnosti ")
#
#plt.plot(time,rho10re, label = 'realni deo')
#plt.plot(time,rho10im, label = 'imaginarni deo')
#plt.xlabel("vreme")
#plt.ylabel("rho10")
#plt.title("Grafik zavisnosti rho_10 od vremena")