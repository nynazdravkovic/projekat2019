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
I = complex(0,1)

matricaGustineIm = []
matricaGustineRe = []
a = []
b = []
c = []
e = []
f = []
a1 = []
t = 0.
h = 0.02

f1 = Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*op*ro02 - 0.5*op*ro20)
f2 = -gamma01*ro01 - 1.0*I*(0.5*oc*ro02 - 0.5*op*ro21 - ro01*(dc - dp))
f3 = -gamma02*ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01 + 0.5*op*ro00 - 0.5*op*ro22)
f4 = -gamma10*ro10 - 1.0*I*(-0.5*oc*ro20 + 0.5*op*ro12 + ro10*(dc - dp))
f5 = -Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21)
f6 = -gamma12*ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + 0.5*op*ro10 + ro12*(dc - dp))
f7 = -gamma20*ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10 - 0.5*op*ro00 + 0.5*op*ro22)
f8 = -gamma21*ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - 0.5*op*ro01 - ro21*(dc - dp))
f9 = -Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21 - 0.5*op*ro02 + 0.5*op*ro20)

 
while (t<4):
    k11 = Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*op*ro02 - 0.5*op*ro20)
    k21 = -gamma01*ro01 - 1.0*I*(0.5*oc*ro02 - 0.5*op*ro21 - ro01*(dc - dp))
    k31 = -gamma02*ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01 + 0.5*op*ro00 - 0.5*op*ro22)
    k41 = -gamma10*ro10 - 1.0*I*(-0.5*oc*ro20 + 0.5*op*ro12 + ro10*(dc - dp)) 
    k51 = -Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21)
    k61 = -gamma12*ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + 0.5*op*ro10 + ro12*(dc - dp))
    k71 = -gamma20*ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10 - 0.5*op*ro00 + 0.5*op*ro22)
    k81 = -gamma21*ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - 0.5*op*ro01 - ro21*(dc - dp))
    k91 = -Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21 - 0.5*op*ro02 + 0.5*op*ro20)
#    ft2 = t + (h/2.)
    ro002 = ro00+(h/2.)*k11
    ro012 = ro01+(h/2.)*k21
    ro022 = ro02+(h/2.)*k31
    ro102 = ro10+(h/2.)*k41
    ro112 = ro11+(h/2.)*k51
    ro122 = ro12+(h/2.)*k61
    ro202 = ro20+(h/2.)*k71
    ro212 = ro21+(h/2.)*k81
    ro222 = ro22+(h/2.)*k91
    k12 = Gamma01*ro112 + Gamma02*ro222 - 1.0*I*(0.5*op*ro022 - 0.5*op*ro202)
    k22 = -gamma01*ro012 - 1.0*I*(0.5*oc*ro022 - 0.5*op*ro212 - ro012*(dc - dp))
    k32 = -gamma02*ro022 - 1.0*I*(dp*ro022 + 0.5*oc*ro012 + 0.5*op*ro002 - 0.5*op*ro222)
    k42 = -gamma10*ro102 - 1.0*I*(-0.5*oc*ro202 + 0.5*op*ro122 + ro102*(dc - dp)) 
    k52 = -Gamma01*ro112 + Gamma12*ro222 - 1.0*I*(0.5*oc*ro122 - 0.5*oc*ro212)
    k62 = -gamma12*ro122 - 1.0*I*(dp*ro122 + 0.5*oc*ro112 - 0.5*oc*ro222 + 0.5*op*ro102 + ro122*(dc - dp))
    k72 = -gamma20*ro202 - 1.0*I*(-dp*ro202 - 0.5*oc*ro102 - 0.5*op*ro002 + 0.5*op*ro222)
    k82 = -gamma21*ro212 - 1.0*I*(-dp*ro212 - 0.5*oc*ro112 + 0.5*oc*ro222 - 0.5*op*ro012 - ro212*(dc - dp))
    k92 = -Gamma02*ro222 - Gamma12*ro222 - 1.0*I*(-0.5*oc*ro122 + 0.5*oc*ro212 - 0.5*op*ro022 + 0.5*op*ro202)
#    ft3 = t + (h/2.)
    ro003 = ro00+(h/2.)*k12
    ro013 = ro01+(h/2.)*k22
    ro023 = ro02+(h/2.)*k32
    ro103 = ro10+(h/2.)*k42
    ro113 = ro11+(h/2.)*k52
    ro123 = ro12+(h/2.)*k62
    ro203 = ro20+(h/2.)*k72
    ro213 = ro21+(h/2.)*k82
    ro223 = ro22+(h/2.)*k92
    k13 = Gamma01*ro113 + Gamma02*ro223 - 1.0*I*(0.5*op*ro023 - 0.5*op*ro203)
    k23 = -gamma01*ro013 - 1.0*I*(0.5*oc*ro023 - 0.5*op*ro213 - ro013*(dc - dp))
    k33 = -gamma02*ro023 - 1.0*I*(dp*ro023 + 0.5*oc*ro013 + 0.5*op*ro003 - 0.5*op*ro223)
    k43 = -gamma10*ro103 - 1.0*I*(-0.5*oc*ro203 + 0.5*op*ro123 + ro103*(dc - dp)) 
    k53 = -Gamma01*ro113 + Gamma12*ro223 - 1.0*I*(0.5*oc*ro123 - 0.5*oc*ro213)
    k63 = -gamma12*ro123 - 1.0*I*(dp*ro123 + 0.5*oc*ro113 - 0.5*oc*ro223 + 0.5*op*ro103 + ro123*(dc - dp))
    k73 = -gamma20*ro203 - 1.0*I*(-dp*ro203 - 0.5*oc*ro103 - 0.5*op*ro003 + 0.5*op*ro223)
    k83 = -gamma21*ro213 - 1.0*I*(-dp*ro213 - 0.5*oc*ro113 + 0.5*oc*ro223 - 0.5*op*ro013 - ro213*(dc - dp))
    k93 = -Gamma02*ro223 - Gamma12*ro223 - 1.0*I*(-0.5*oc*ro123 + 0.5*oc*ro213 - 0.5*op*ro023 + 0.5*op*ro203)
#    ft4 = t + (h/2.)
    ro004 = ro00+(h/2.)*k13
    ro014 = ro01+(h/2.)*k23
    ro024 = ro02+(h/2.)*k33
    ro104 = ro10+(h/2.)*k43
    ro114 = ro11+(h/2.)*k53
    ro124 = ro12+(h/2.)*k63
    ro204 = ro20+(h/2.)*k73
    ro214 = ro21+(h/2.)*k83
    ro224 = ro22+(h/2.)*k93
    k14 = Gamma01*ro114 + Gamma02*ro224 - 1.0*I*(0.5*op*ro024 - 0.5*op*ro204)
    k24 = -gamma01*ro014 - 1.0*I*(0.5*oc*ro024 - 0.5*op*ro214 - ro014*(dc - dp))
    k34 = -gamma02*ro024 - 1.0*I*(dp*ro024 + 0.5*oc*ro014 + 0.5*op*ro004 - 0.5*op*ro224)
    k44 = -gamma10*ro104 - 1.0*I*(-0.5*oc*ro204 + 0.5*op*ro124 + ro104*(dc - dp)) 
    k54 = -Gamma01*ro114 + Gamma12*ro224 - 1.0*I*(0.5*oc*ro124 - 0.5*oc*ro214)
    k64 = -gamma12*ro124 - 1.0*I*(dp*ro124 + 0.5*oc*ro114 - 0.5*oc*ro224 + 0.5*op*ro104 + ro124*(dc - dp))
    k74 = -gamma20*ro204 - 1.0*I*(-dp*ro204 - 0.5*oc*ro104 - 0.5*op*ro004 + 0.5*op*ro224)
    k84 = -gamma21*ro214 - 1.0*I*(-dp*ro214 - 0.5*oc*ro114 + 0.5*oc*ro224 - 0.5*op*ro014 - ro214*(dc - dp))
    k94 = -Gamma02*ro224 - Gamma12*ro224 - 1.0*I*(-0.5*oc*ro124 + 0.5*oc*ro214 - 0.5*op*ro024 + 0.5*op*ro204)
    t = t + (h/2.)
    ro00 = ro00 + (h/6.)*(k11+(2.*k12)+(2.*k13)+k14)
    ro01 = ro01 + (h/6.)*(k21+(2.*k22)+(2.*k23)+k24)
    ro02 = ro02 + (h/6.)*(k31+(2.*k32)+(2.*k33)+k34)
    ro10 = ro10 + (h/6.)*(k41+(2.*k42)+(2.*k43)+k44)
    ro11 = ro11 + (h/6.)*(k51+(2.*k52)+(2.*k53)+k54)
    ro12 = ro12 + (h/6.)*(k61+(2.*k62)+(2.*k63)+k64)
    ro20 = ro20 + (h/6.)*(k71+(2.*k72)+(2.*k73)+k74)
    ro21 = ro21 + (h/6.)*(k81+(2.*k82)+(2.*k83)+k84)
    ro22 = ro22 + (h/6.)*(k91+(2.*k92)+(2.*k93)+k94)
    a.append(ro00.real)
    e.append(ro11.real)
    f.append(ro22.real)

    c.append(ro12.imag)
    a1.append(ro12.real)
    b.append(t)

#plotovanje 2d
plt.plot(b,a)
plt.title("Realni deo")
plt.xlabel("t")
plt.ylabel("rho00")
plt.savefig('ro00re.png')
plt.show()
plt.plot(b,e)
plt.title("Realni deo")
plt.xlabel("t")
plt.ylabel("rho11")
plt.savefig('ro00re.png')
plt.show()
plt.plot(b,f)
plt.title("Realni deo")
plt.xlabel("t")
plt.ylabel("rho22")
plt.savefig('ro00re.png')
plt.show()
 
plt.plot(b,a)
plt.title("Realni deo")
plt.xlabel("t")
plt.ylabel("rho01")
plt.savefig('ro01re.png')
plt.show()
   
plt.plot(b,c)
plt.title("Imaginarni deo")
plt.xlabel("t")
plt.ylabel("rho01")
plt.savefig('ro01im.png')

plt.show()