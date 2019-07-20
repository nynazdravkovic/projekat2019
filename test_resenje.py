# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 17:30:26 2019

@author: nina
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#u ovom kodu plotujem samo jednacine
o = 0.5
j = complex(0,1)
ro01_im=[]
ro00_im=[]
ro10_im=[]
ro11_im=[]
ro11_re=[]
ro00_re=[]
ro01_re=[]
ro10_re=[]
delta=np.linspace(-500,500,1000)
t=100

for i in range (1000):
    d = (i-500)/500    
    if d == 0:
        rho01=0
        rho11=0
        rho00=0
        rho10=0
    else:    
        rho01= -o/(2*d) - o*np.exp(-j*t*np.sqrt(4*d**2 + 2*o**2)/2)/(2*d - np.sqrt(2)*np.sqrt(2*d**2 + o**2)) - o*np.exp(j*t*np.sqrt(4*d**2 + 2*o**2)/2)/(2*d + np.sqrt(2)*np.sqrt(2*d**2 + o**2)) 
        rho10= -o/(2*d) - o*np.exp(-j*t*np.sqrt(4*d**2 + 2*o**2)/2)/(2*d + np.sqrt(2)*np.sqrt(2*d**2 + o**2)) - o*np.exp(j*t*np.sqrt(4*d**2 + 2*o**2)/2)/(2*d - np.sqrt(2)*np.sqrt(2*d**2 + o**2)) 
        rho11= 1 + np.exp(-j*t*np.sqrt(4*d**2 + 2*o**2)/2) + np.exp(j*t*np.sqrt(4*d**2 + 2*o**2)/2)
        rho00= -np.exp(-j*t*np.sqrt(4*d**2 + 2*o**2)/2) + np.exp(j*t*np.sqrt(4*d**2 + 2*o**2)/2)
    ro01_im.append(rho01.imag)
    ro01_re.append(rho01.real)
    ro10_im.append(rho10.imag)
    ro10_re.append(rho10.real)
        
plt.plot(delta,ro01_im)
plt.title("Imaginarni deo")
plt.xlabel("delta")
plt.ylabel("rho01")
plt.show()
plt.plot(delta,ro10_im)
plt.title("Imaginarni deo")
plt.xlabel("delta")
plt.ylabel("rho10")
plt.show()
plt.plot(delta,ro01_re)
plt.title("Realni deo")
plt.xlabel("delta")
plt.ylabel("rho01")
plt.show()
plt.plot(delta,ro10_re)
plt.title("Realni deo")
plt.xlabel("delta")
plt.ylabel("rho10")
plt.show()
