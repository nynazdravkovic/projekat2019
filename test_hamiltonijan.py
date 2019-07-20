# -*- coding: utf-8 -*-
"""
Created on Fri Jun 28 14:31:37 2019

@author: nina
"""

import sympy
import numpy
from sympy import *
import math
#dobijanje hamiltonijana

w=sympy.Symbol('w')
o=sympy.Symbol('o')
h=sympy.Symbol('h')
fi=sympy.Symbol('fi')
E1=sympy.Symbol('E1')
E0=sympy.Symbol('E0')
t=sympy.Symbol('t')
d=sympy.Symbol('d')
I = complex(0,1)

uc = [[1,0],[0,sympy.exp(I*(w*t + fi))]]
u = [[1,0],[0,sympy.exp(-I*(w*t + fi))]]
ud = [[0,0],[0,I*w*sympy.exp(I*(w*t + fi))]]
H = [[E0, -h/2*o*sympy.exp(I*(w*t + fi))],[-h/2*o*sympy.exp(-I*(w*t + fi)),E1]]

result0 = [[0,0],[0,0]]


def ucH():
    for i in range(len(uc)):
       for j in range(len(H[0])):
           for k in range(len(H)):
               result0[i][j] += uc[i][k] *H[k][j]
    b = result0
    print(b)
    return b

ucH=ucH()
result1 = [[0,0],
         [0,0]]
def ucHu():
    for i in range(len(ucH)):
       for j in range(len(u[0])):
           for k in range(len(u)):
               result1[i][j] += ucH[i][k] *u[k][j]
    a = result1
    print(a)
    return a
ucHu = ucHu()
result2 = [[0,0],
         [0,0]]
def udu():
    for i in range(len(ud)):
       for j in range(len(u[0])):
           for k in range(len(ud)):
               result2[i][j] += I*h*ud[i][k] *u[k][j]
    c = result2
    print(c)
    return c
udu=udu()
rezultat=[[0,0],
         [0,0]]
for i in range(len(ucHu)):
   for j in range(len(ucHu[0])):
       rezultat[i][j] = ucHu[i][j] + udu[i][j]

for r in rezultat:
    print(r)