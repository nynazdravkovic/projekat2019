# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:23:48 2019

@author: nina
"""
import matplotlib.pyplot as plt
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
dp = 0.1
dc = 0.1
I = complex(0,1)
ro = [ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22]
k = np.full((4,9),(complex(0,0)))

f1 = Gamma01*ro[4] + Gamma02*ro[8] - 1.0*I*(0.5*op*ro[2] - 0.5*op*ro[6])
f2 = -gamma01*ro[1] - 1.0*I*(0.5*oc*ro[2] - 0.5*op*ro[7] - ro[1]*(dc - dp))
f3 = -gamma02*ro[2] - 1.0*I*(dp*ro[2] + 0.5*oc*ro[1] + 0.5*op*ro[0] - 0.5*op*ro[8])
f4 = -gamma10*ro[3] - 1.0*I*(-0.5*oc*ro[6] + 0.5*op*ro[5] + ro[3]*(dc - dp))
f5 = -Gamma01*ro[4] + Gamma12*ro[8] - 1.0*I*(0.5*oc*ro[5] - 0.5*oc*ro[7])
f6 = -gamma12*ro[5] - 1.0*I*(dp*ro[5] + 0.5*oc*ro[4] - 0.5*oc*ro[8] + 0.5*op*ro[3] + ro[5]*(dc - dp))
f7 = -gamma20*ro[6] - 1.0*I*(-dp*ro[6] - 0.5*oc*ro[3] - 0.5*op*ro[0] + 0.5*op*ro[8])
f8 = -gamma21*ro[7] - 1.0*I*(-dp*ro[7] - 0.5*oc*ro[4] + 0.5*oc*ro[8] - 0.5*op*ro[1] - ro[7]*(dc - dp))
f9 = -Gamma02*ro[8] - Gamma12*ro[8] - 1.0*I*(-0.5*oc*ro[5] + 0.5*oc*ro[7] - 0.5*op*ro[2] + 0.5*op*ro[6])
c= []
F = [f1,f2,f3,f4,f5,f6,f7,f8,f9]
t = 0.
brojPonavljanja = 100
h = 0.0005
print(F[1])
matricaGustineIm = np.zeros((9,brojPonavljanja))
matricaGustineRe = np.zeros((9,brojPonavljanja))
for x in range(brojPonavljanja):
    for j in range(4):
        print('/')
        for i in range(9):
            k[j][i] = F[i]
            F[i] = ro[i]+h*0.5*k[j][i]
        print(F)
    for i in range(9):
        F[i] += (1/6)*(k[0][i]+2*k[1][i]+2*k[2][i]+k[3][i])
    for i in range(9):    
         matricaGustineIm[i][x] = F[i].imag
         matricaGustineRe[i][x] = F[i].real
    c.append(t)
    t+=h/2
for i in range(9):
    plt.plot(c,matricaGustineRe[i])
    plt.title("Realni deo")
    plt.xlabel("t")
    plt.show()
for i in range(9):
    plt.plot(c,matricaGustineIm[i])
    plt.title("Imaginarni deo")
    plt.xlabel("t")
    plt.show()
#plt.title("Realni deo")
#plt.xlabel("t")
#plt.ylabel("rho20")
#plt.savefig('ro01re_3n.png')
#plt.show()
#   
#plt.plot(c,a)
#plt.title("Imaginarni deo")
#plt.xlabel("t")
#plt.ylabel("rho20")
#plt.savefig('ro01im_3n.png')
#plt.show()
#
#plt.plot(c,e)
#plt.title("Realni deo")
#plt.xlabel("t")
#plt.ylabel("rho00")
#plt.savefig('ro00re_3n.png')
#plt.show()
