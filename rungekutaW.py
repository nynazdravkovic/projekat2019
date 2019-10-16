# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 00:14:14 2019

@author: nina
"""

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
gamma30 = 1.
gamma31 = 1.
gamma32 = 1.
gamma33 = 1.
op = 0.5
#oc = 0.5
dp = 0.
#dc = 0.
oc1 = 0.5
oc2 = 0.5
#oc = 0.5
dc1 = 0.
dc2 = 0.
#dc = 0.
I = complex(0,1)

matricaGustineIm = []
matricaGustineRe = []
rho00re = []
rho01re = []
rho02re = []
rho03re = []
rho10re = []
rho11re = []
rho12re = []
rho13re = []
rho20re = []
rho21re = []
rho22re = []
rho23re = []
rho30re = []
rho31re = []
rho32re = []
rho33re = []

rho00im = []
rho01im = []
rho02im = []
rho03im = []
rho10im = []
rho11im = []
rho12im = []
rho13im = []
rho20im = []
rho21im = []
rho22im = []
rho23im = []
rho30im = []
rho31im = []
rho32im = []
rho33im = []
t = 0.
h = 0.02
ro00re = []
time = []
nultiStepen=[[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1],
                 [-I*oc1,-gamma01+I*dc1,0,0,0,I*oc1,0,0,I*oc2,0,0,0,0,0,0,0],
                 [-I*oc2,0,-gamma02+I*dc2,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0,0],
                 [0,0,0,I*dp-gamma03,0,0,0,I*oc1,0,0,0,I*oc2,0,0,0,0],
                 [I*oc1,0,0,0,-gamma10-I*dc1,-I*oc1,-oc2*I,0,0,0,0,0,0,0,0,0],
                 [0,I*oc1,0,0,-I*oc1,-Gamma01,0,0,0,0,Gamma12,Gamma13,0,0,0,0],
                 [0,0,I*oc1,0,-I*oc2,0,dc2-dc1-gamma12,0,0,0,0,0,0,0,0,0],
                 [0,0,0,I*oc1,0,0,0,dp-dc1-gamma13,0,0,0,0,0,0,0,0],
                 [-I*oc2,0,0,0,0,0,0,0,-gamma20-I*dc2,-I*oc1,-I*oc2,0,0,0,0,0],
                 [0,I*oc2,0,0,0,0,0,0,-I*oc1,-gamma21+I*(dc1-dc2),0,0,0,0,0,0],
                 [0,0,I*oc2,0,0,0,0,0,-I*oc2,0,-Gamma02-Gamma12,0,0,0,0,Gamma23],
                 [0,0,0,oc2*I,0,0,0,0,0,0,0,-gamma23+I*(dp-dc2),0,0,0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-gamma30-I*dp,-I*oc1,-I*oc2,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc1,-gamma31-I*(dp-dc1),0,0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,-I*oc2,0,-gamma32-I*(dp-dc2),0],
                 [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,Gamma03]]

nule = np.zeros(16)
nule[0] = 1
pocetniUslovi = np.linalg.solve(nultiStepen,nule)
ro00 = pocetniUslovi[0]
ro01 = pocetniUslovi[1]
ro02 = pocetniUslovi[2]
ro03 = pocetniUslovi[3]
ro10 = pocetniUslovi[4]
ro11 = pocetniUslovi[5]
ro12 = pocetniUslovi[6]
ro13 = pocetniUslovi[7]
ro20 = pocetniUslovi[8]
ro21 = pocetniUslovi[9]
ro22 = pocetniUslovi[10]
ro23 = pocetniUslovi[11]
ro30 = pocetniUslovi[12]
ro31 = pocetniUslovi[13]
ro32 = pocetniUslovi[14]
ro33 = pocetniUslovi[15]
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
    f1 = Gamma01*ro[5] + Gamma02*ro[10] + Gamma03*ro[15] - 1.0*I*(oc1*ro[1] - oc1*ro[4] + oc2*ro[2] - oc2*ro[8] + op*ro[3] - op*ro[12])
    f2=-gamma01*ro[1] - 1.0*I*(-dc1*ro[1] + oc1*ro[0] - oc1*ro[5] - oc2*ro[11] - op*ro[13])
    f3=-gamma02*ro[2] - 1.0*I*(-dc2*ro[2] - oc1*ro[6] + oc2*ro[0] - oc2*ro[10] - op*ro[14])
    f4=-gamma03*ro[3] - 1.0*I*(-dp*ro[3] - oc1*ro[7] - oc2*ro[11] + op*ro[0] - op*ro[15])
    f5=-gamma10*ro[4] - 1.0*I*(dc1*ro[4] - oc1*ro[0] + oc1*ro[5] + oc2*ro[6] + op*ro[7])
    f6=-Gamma01*ro[5] + Gamma12*ro[10] + Gamma13*ro[15] - 1.0*I*(-oc1*ro[1] + oc1*ro[4])
    f7=-gamma12*ro[6] - 1.0*I*(dc1*ro[6] - dc2*ro[6] - oc1*ro[2] + oc2*ro[4])
    f8=-gamma13*ro[7] - 1.0*I*(dc1*ro[7]- dp*ro[7] - oc1*ro[12] + op*ro[4])
    f9=-gamma20*ro[8] - 1.0*I*(dc2*ro[8] + oc1*ro[9] - oc2*ro[0] + oc2*ro[10] + op*ro[11])
    f10=-gamma21*ro[9] - 1.0*I*(-dc1*ro[9] + dc2*ro[9] + oc1*ro[8] - oc2*ro[1])
    f11=-Gamma02*ro[10] - Gamma12*ro[10] + Gamma23*ro[15] - 1.0*I*(-oc2*ro[2] + oc2*ro[8])
    f12=-gamma23*ro[11] - 1.0*I*(dc2*ro[11] - dp*ro[11] - oc2*ro[3] + op*ro[8])
    f13=-gamma30*ro[12] - 1.0*I*(dp*ro[12] + oc1*ro[13] + oc2*ro[14] - op*ro[0] + op*ro[15])
    f14=-gamma31*ro[13] - 1.0*I*(-dc1*ro[13] + dp*ro[13] + oc1*ro30 - op*ro01)
    f15=-gamma32*ro[14] - 1.0*I*(-dc2*ro[14] + dp*ro[14] + oc2*ro[12] - op*ro[2])
    f16=-Gamma02*ro[15] - Gamma12*ro[15] - 1.0*I*(-op*ro[3] + op*ro[12])
    return([f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15,f16])

while (t<4):
    K1=np.array(mnozenjeK(pocetniUslovi))*h
    K2=np.array(mnozenjeK(pocetniUslovi+K1/2))*h
    K3=np.array(mnozenjeK(pocetniUslovi+K2/2))*h
    K4=np.array(mnozenjeK(pocetniUslovi+K3))*h
    t = t + (h)  
    for j in range(16):
        pocetniUslovi[j] += (1/6.)*(K1[j]+(2.*K2[j])+(2.*K3[j])+K4[j])
    time.append(t)  
    rho00re.append(pocetniUslovi[0].real)
    rho01re.append(pocetniUslovi[1].real)
    rho02re.append(pocetniUslovi[2].real)
    rho03re.append(pocetniUslovi[3].real)
    rho10re.append(pocetniUslovi[4].real)
    rho11re.append(pocetniUslovi[5].real)
    rho12re.append(pocetniUslovi[6].real)
    rho13re.append(pocetniUslovi[7].real)
    rho20re.append(pocetniUslovi[8].real)
    rho21re.append(pocetniUslovi[9].real)
    rho22re.append(pocetniUslovi[10].real)
    rho23re.append(pocetniUslovi[11].real)
    rho30re.append(pocetniUslovi[12].real)
    rho31re.append(pocetniUslovi[13].real)
    rho32re.append(pocetniUslovi[14].real)
    rho33re.append(pocetniUslovi[15].real)
    
    rho00im.append(pocetniUslovi[0].imag)
    rho01im.append(pocetniUslovi[1].imag)
    rho02im.append(pocetniUslovi[2].imag)
    rho03im.append(pocetniUslovi[3].imag)
    rho10im.append(pocetniUslovi[4].imag)
    rho11im.append(pocetniUslovi[5].imag)
    rho12im.append(pocetniUslovi[6].imag)
    rho13im.append(pocetniUslovi[7].imag)
    rho20im.append(pocetniUslovi[8].imag)
    rho21im.append(pocetniUslovi[9].imag)
    rho22im.append(pocetniUslovi[10].imag)
    rho23im.append(pocetniUslovi[11].imag)
    rho30im.append(pocetniUslovi[12].imag)
    rho31im.append(pocetniUslovi[13].imag)
    rho32im.append(pocetniUslovi[14].imag)
    rho33im.append(pocetniUslovi[15].imag)
    
plt.plot(time,rho00re, 'g:', label = 'rho00')
#plt.plot(time,rho01re,'b:' , label = 'rho01')
#plt.plot(time,rho02re, 'r:', label = 'rho02')
#plt.plot(time,rho10re,'m:' , label = 'rho10')
plt.plot(time,rho11re, 'g-', label = 'rho11')
#plt.plot(time,rho12re, 'b-', label = 'rho12')
#plt.plot(time,rho20re, 'r-', label = 'rho20')
#plt.plot(time,rho21re, 'm', label = 'rho21')
plt.plot(time,rho22re, 'm--', label = 'rho22')
plt.plot(time,rho33re, 'm:', label = 'rho33')


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
a = []
for i in range(len(time)):
    a.append(rho11re[i]+rho22re[i]+rho33re[i]+rho00re[i])
plt.plot(time,a, label = 'realni deo')
#plt.legend()
#plt.xlabel("vreme")
#plt.ylabel("rho02")
plt.show()
#plt.title("Grafik zavisnosti ")
#
#plt.plot(time,rho10re, label = 'realni deo')
#plt.plot(time,rho10im, label = 'imaginarni deo')
#plt.xlabel("vreme")
#plt.ylabel("rho10")
#plt.title("Grafik zavisnosti rho_10 od vremena")