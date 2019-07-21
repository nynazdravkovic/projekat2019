# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 17:20:07 2019

@author: nina
"""
import sympy
import numpy
from sympy import *
import math
import scipy.sparse.linalg
import matplotlib.pyplot as plt


#ovde izvodimo matricu gustine za sistem sa dva nivoa i dva lasera

w=sympy.Symbol('w')
h=sympy.Symbol('h')
#d=sympy.Symbol('d')
o=sympy.Symbol('o')
I = sympy.I
ro00=sympy.Symbol('ro00')
ro01=sympy.Symbol('ro01')
ro10=sympy.Symbol('ro10')
ro11=sympy.Symbol('ro11')
gamma00=1
gamma01=1
gamma10=1
gamma11=1

brojac=[]
imaginarnaGustina1=[]
realnaGustina1=[]
imaginarnaGustina2=[]
realnaGustina2=[]


ro=[[ro00, ro01],[ro10, ro11]]
gamma=[[gamma00, gamma01],[gamma10, gamma11]]
Gamma=[[1, 1, 1],[ 1,1,1],[1,1,1]]


result0 = [[0,0],[0,0]]
result1 = [[0,0],[0,0]]
result = [[0,0],[0,0]]

d = sympy.Symbol('d')
H = [[0, -h*o/2],[-h*o/2, h*d]]
for i in range(len(H)):
   for j in range(len(ro[0])):
       for k in range(len(ro)):
           result0[i][j] += H[i][k] * ro[k][j]

for i in range(len(ro)):
   for j in range(len(H[0])):
       for k in range(len(H)):
           result1[i][j] += ro[i][k] * H[k][j]
           
for i in range (2):
    for j in range (2):
        result[i][j]=-I*(result0[i][j]-result1[i][j])/h
        if i!=j:
            result[i][j]+=-gamma[i][j]*ro[i][j]
        else:
            if i == 0:
                result[i][j]=result[i][j]+Gamma[0][1]*ro[1][1]
            else:
                result[i][j]=result[i][j]-Gamma[0][1]*ro[1][1]
            
#    print(result)    
      
for i in range (2):
    for j in range (2):
        locals()["jna"+str(i)+str(j)] = result[i][j]
#            print(result[i][j])
jna = ro00+ro11-1
for i in range (2):
        for j in range (2):
            locals()["jna"+str(i)+str(j)] = result[i][j]
            print(result[i][j])
#jednacina = sympy.solve([jna01,jna10,jna11,jna],[ro00,ro01,ro10,ro11])
#print(jednacina)