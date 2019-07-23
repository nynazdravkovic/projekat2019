# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 17:49:46 2019

@author: nina
"""
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

n = 100
deltax = 0.01
ro12im = np.loadtxt("ro12im.csv")
ro12re = np.loadtxt("ro12re.csv")
dEdt = np.empty((n,n),dtype = 'complex')
matrica = np.zeros((n,n))
#A mi je (f[i+1]-f[i-1])
A = np.zeros(n)
E = np.linspace(0,1,100)
z0 = np.linspace(0,n,n)
for i in range (n):
    for j in range(n):
        if(i==j):
            if(i==0):
                matrica[0][n-1] = 1
            else:
                matrica[i-1][j] = 1

for i in range (n):
    for j in range(n):
        if(i==j):
            if(i==0):
                matrica[0][n-1] = -1
            else:
                matrica[i][j-1] = -1

for i in range(n):
    for j in range(n):
        A[i]+=matrica[i][j]*E[j]
def sistem(t, z):
    for i in range(n):
        E[i]=z[i]
        dEdt[i] = -1*A[i]*ro12im[i]*0.5/deltax
    return(dEdt)
t=np.linspace(0,5,100)
z = solve_ivp(sistem, (0, 5), z0)

