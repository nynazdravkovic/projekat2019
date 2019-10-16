# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:06 2019

@author: nina
"""
#

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.integrate 
from scipy.integrate import quad, dblquad, tplquad, nquad
import scipy.special

#I=quad(functionF, opseg1FLOAt, opseg2FLOAT, args=(rho,phi,z))
#I=tplquad(functionF, opseg1F, opseg2F, lowY, upY, lowZ, upZ)

#DEFINISANJE
m0=0.665*9.1*10**(-31)
me=0.38*9.1*10**(-31)
hb=6.62*10**(-34)
R=10*10**(-9)
L=10*10**(-9)

q=1.6*10**(-19)
#rho=1
#fi=1
#z=1

#rho=np.linspace(0,10,10)
#fi=np.linspace(0,10,10)
#z=np.linspace(0,10,10)

omega0rho=hb/(me*R**2)
omega0z=4*hb/(me*L**2)

krho=math.sqrt(m0*omega0rho/hb)
kz=math.sqrt(m0*omega0z/hb)

def psi0(rho, fi, z):
    preE=math.sqrt(kz*krho**2)/math.pi**(3/2)
    posleE=math.e**(-1/2*(kz**2*z**2+krho**2*rho**2))
    return preE*posleE
def psi1(rho, fi, z):
    preE=math.sqrt(kz*krho**2)/(2*math.sqrt(2)*math.pi**(3/2))*(4*kz**2*z**2-2)
    posleE=math.e**(-1/2*(kz**2*z**2+krho**2*rho**2))
    return preE*posleE
def psi2(rho, fi, z):
    preE=math.sqrt(kz*krho**2)/(4*math.sqrt(3)*math.pi**(3/2))*(8*kz*z**3-12*kz*z)
    posleE=math.e**(-1/2*(kz**2*z**2+krho**2*rho**2))
    return preE*posleE
def psi2ct(rho, fi, z):
    return np.transpose((np.conj(psi2(rho, fi, z))))
def psi2ctpsi0(rho, fi, z):
    return(psi2ct(rho, fi, z)*psi0(rho, fi, z))
    

#PROVERITI LAMBDA I GRANICNE VREDNOSTI I DA LI IMA NECEGA IZMEDJU OVE 2 FUNKCIJE 
#prvo su za prvi integral geanicne, drugo za drugi, trece za treci
#mi02=q*tplquad(psi2ctpsi0, 0, 10,  lambda z: 0, lambda z:10, lambda z, rho: 0, lambda z, rho: 10)


#mi02=nquad(psi2ctpsi0, [(-math.sqrt(2),math.sqrt(2)), (-math.pi/4,math.pi/4), (0,L)])
mi02=nquad(psi2ctpsi0, [(0,np.inf), (-math.pi,math.pi), (-np.inf, np.inf)])

print(q*mi02[0])

mii= -math.e*kz**2*L**3/(4*math.sqrt(3)*math.pi**2)*(math.e**(krho**2*R**2)-1)*math.e**(-1/4*(4*kz**2*L**2+krho**2*R**2))   
#mii=math.sqrt(kz*krho**2)/(2*math.sqrt(2)*math.pi**(3/2))*(4*krho**2*rho**2*math.cos(fi)-2)*math.e**(-1/2*(kz**2*z**2+krho**2*rho**2))   
print(mii)
