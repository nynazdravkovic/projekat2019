# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:06 2019

@author: nina
"""
#
#W konfiguracija
import sympy
import numpy as np
from sympy import *
from sympy.solvers.solvers import solve_linear_system
import scipy.sparse.linalg
import scipy.sparse
#resenje u brojevima
#dc=0
#op=0.1
#oc=1
##resenje u opstim brojevima
dc=sympy.Symbol('dc')
dp=sympy.Symbol('dp')
oc=sympy.Symbol('oc')
op=sympy.Symbol('op')

ro00=sympy.Symbol('ro00')
ro01=sympy.Symbol('ro01')
ro02=sympy.Symbol('ro02')
ro03=sympy.Symbol('ro03')
ro10=sympy.Symbol('ro10')
ro11=sympy.Symbol('ro11')
ro12=sympy.Symbol('ro12')
ro13=sympy.Symbol('ro13')
ro20=sympy.Symbol('ro20')
ro21=sympy.Symbol('ro21')
ro22=sympy.Symbol('ro22')
ro23=sympy.Symbol('ro23')
ro30=sympy.Symbol('ro30')
ro31=sympy.Symbol('ro31')
ro32=sympy.Symbol('ro32')
ro33=sympy.Symbol('ro33')

gamma00=sympy.Symbol('gamma00')
gamma01=sympy.Symbol('gamma01')
gamma02=sympy.Symbol('gamma02')
gamma03=sympy.Symbol('gamma03')
gamma10=sympy.Symbol('gamma10')
gamma11=sympy.Symbol('gamma11')
gamma12=sympy.Symbol('gamma12')
gamma13=sympy.Symbol('gamma13')
gamma20=sympy.Symbol('gamma20')
gamma21=sympy.Symbol('gamma21')
gamma22=sympy.Symbol('gamma22')
gamma23=sympy.Symbol('gamma23')
gamma30=sympy.Symbol('gamma30')
gamma31=sympy.Symbol('gamma31')
gamma32=sympy.Symbol('gamma32')
gamma33=sympy.Symbol('gamma33')

Gamma00=sympy.Symbol('Gamma00')
Gamma01=sympy.Symbol('Gamma01')
Gamma02=sympy.Symbol('Gamma02')
Gamma03=sympy.Symbol('Gamma03')
Gamma10=sympy.Symbol('Gamma10')
Gamma11=sympy.Symbol('Gamma11')
Gamma12=sympy.Symbol('Gamma12')
Gamma13=sympy.Symbol('Gamma13')
Gamma20=sympy.Symbol('Gamma20')
Gamma21=sympy.Symbol('Gamma21')
Gamma22=sympy.Symbol('Gamma22')
Gamma23=sympy.Symbol('Gamma23')
Gamma30=sympy.Symbol('Gamma30')
Gamma31=sympy.Symbol('Gamma31')
Gamma32=sympy.Symbol('Gamma32')
Gamma33=sympy.Symbol('Gamma33')
h = sympy.Symbol('h')
ro=[[ro00, ro01, ro02,ro03],[ro10, ro11, ro12,ro13],[ro20, ro21, ro22,ro23],[ro30,ro31,ro32,ro33]]
gamma=[[gamma00, gamma01, gamma02,gamma03],[gamma10, gamma11, gamma12,gamma13],[gamma20, gamma21, gamma22,gamma23],[gamma30,gamma31,gamma32,gamma33]]
Gamma=[[Gamma00, Gamma01, Gamma02,Gamma03],[Gamma10, Gamma11, Gamma12,Gamma13],[Gamma20, Gamma21, Gamma22,Gamma23],[Gamma30,Gamma31,Gamma32,Gamma33]]
ii = sympy.I

result0 = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
result1 = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
result = [[0,0,0,0],
         [0,0,0,0],
         [0,0,0,0],
         [0,0,0,0]]
brojac=[]
oc1 = sympy.Symbol('oc1')
oc2 = sympy.Symbol('oc2')
op = sympy.Symbol('op')
dc2 = sympy.Symbol('dc2')
dp = sympy.Symbol('dp')
dc1 = sympy.Symbol('dc1')
H=[[0,-oc1,-oc2,-op],[-oc1,dc1,0,0],[-oc2,0,dc2,0],[-op,0,0,dp]]

fi = 0
for i in range(len(H)):
   for j in range(len(ro[0])):
       for k in range(len(ro)):
           result0[i][j] += H[i][k] * ro[k][j]

for i in range(len(ro)):
   for j in range(len(H[0])):
       for k in range(len(H)):
           result1[i][j] += ro[i][k] * H[k][j]
           
for i in range (4):
    for j in range (4):
        result[i][j]=-complex(0,1)*(result0[i][j]-result1[i][j])
        if i!=j:
            result[i][j]+=-gamma[i][j]*ro[i][j]
        else:
            if i == 0:
                result[i][j]=result[i][j]+Gamma[0][1]*ro[1][1]+Gamma[0][2]*ro[2][2]+Gamma[0][3]*ro[3][3]
            elif i == 1:
                result[i][j]=result[i][j]-Gamma[0][1]*ro[1][1]+Gamma[1][2]*ro[2][2]+Gamma[1][3]*ro[3][3]
            elif i ==2:
                result[i][j]=result[i][j]-Gamma[0][2]*ro[2][2]-Gamma[1][2]*ro[2][2]+Gamma[2][3]*ro[3][3]
            else:
                result[i][j]=result[i][j]-Gamma[0][2]*ro[3][3]-Gamma[1][2]*ro[3][3]
                
                
#print(result)                
for i in range (4):
    for j in range (4):
        locals()["jna"+str(i)+str(j)] = result[i][j]
        print(result[i][j])                

jna = ro00 + ro11 + ro22 - 1
