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
ro00 = 1
ro01 = 0
ro10 = 0
ro11 = 0
o = 0.5
delte = []
d = 0
I = complex(0,1)
matricaGustineIm = []
matricaGustineRe = []
#for j in range(100):
#    a = []
#    c = []
#    b = []
t = 0.
h = 0.02
#    d = (j - 50)/50
#    delte.append(d)
while (t<1):
    
    k11 = ro11 - I*(o*ro01/2. - o*ro10/2.)
    k21 = -ro01 - I*(-d*ro01 + o*ro00/2. - o*ro11/2.)
    k31 = -ro10 - I*(d*ro10 - o*ro00/2. + o*ro11/2.)
    k41 = -ro11 - I*(-o*ro01/2. + o*ro10/2.)
    ft2 = t + (h/2.)
    ro002 = ro00+(h/2.)*k11
    ro012 = ro01+(h/2.)*k21
    ro102 = ro10+(h/2.)*k31
    ro112 = ro11+(h/2.)*k41
    k12 = ro112 - I*(o*ro012/2. - o*ro102/2.)
    k22 = -ro012 - I*(-d*ro012 + o*ro002/2. - o*ro112/2.)
    k32 = -ro102 - I*(d*ro102 - o*ro002/2. + o*ro112/2.)
    k42 = -ro112 - I*(-o*ro012/2. + o*ro102/2.)
    ft3 = t + (h/2.)
    ro003 = ro002+(h/2.)*k12
    ro013 = ro012+(h/2.)*k22
    ro103 = ro102+(h/2.)*k32
    ro113 = ro112+(h/2.)*k42
    k13 = ro113 - I*(o*ro013/2. - o*ro103/2.)
    k23 = -ro013 - I*(-d*ro013 + o*ro003/2. - o*ro113/2.)
    k33 = -ro103 - I*(d*ro103 - o*ro003/2. + o*ro113/2.)
    k43 = -ro113- I*(-o*ro013/2. + o*ro103/2.)
    ft4 = t + (h/2.)
    ro004 = ro003+(h/2.)*k13
    ro014 = ro013+(h/2.)*k23
    ro104 = ro103+(h/2.)*k33
    ro114= ro113+(h/2.)*k43
    k14 = ro114 - I*(o*ro014/2. - o*ro104/2.)
    k24= -ro014- I*(-d*ro014 + o*ro004/2. - o*ro114/2.)
    k34 = -ro104 - I*(d*ro104 - o*ro004/2. + o*ro114/2.)
    k44 = -ro114- I*(-o*ro014/2. + o*ro104/2.)
    t = t + (h/2.)
    ro00 = ro00 + (h/6.)*(k11+(2.*k12)+(2.*k13)+k14)
    ro01 = ro01 + (h/6.)*(k21+(2.*k22)+(2.*k23)+k24)
    ro10 = ro10 + (h/6.)*(k31+(2.*k32)+(2.*k33)+k34)
    ro11 = ro11 + (h/6.)*(k41+(2.*k42)+(2.*k43)+k44)
    a.append(ro10.real)
    c.append(ro10.imag)
    b.append(t)
#matricaGustineRe.append(a)
#matricaGustineIm.append(c)
#print(len(b[1:]))
#print(delte[1:])
#print(len(matricaGustineRe[1:]))

#print(b)


#plotovanje 2d
plt.plot(b,a)
plt.title("Realni deo")
plt.xlabel("t")
plt.ylabel("rho10")
plt.show()
    
#pokusaj plotovanja 3d:
#fig = plt.figure()
#fig.suptitle('Realni deo', fontsize=16)
#ax = fig.add_subplot(111, projection='3d')
#ax.contour3D(b, delte, matricaGustineRe, 50, cmap='binary')
#ax.set_xlabel('t')
#ax.set_ylabel('delta')
#ax.set_zlabel('ro10')


