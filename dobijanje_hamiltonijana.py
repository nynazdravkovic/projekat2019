# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 16:46:42 2019

@author: nina
"""
#OVO SIGURNO RADI I JEDINI JE POTREBAN !!! 24.6.2019.
import sympy
import numpy
from sympy import *
import math

dp=sympy.Symbol('dp')
dc=sympy.Symbol('dc')
fic1=sympy.Symbol('fic1')
fic2=sympy.Symbol('fic2')
fip=sympy.Symbol('fip')
sc1=sympy.Symbol('sc1')
sc2=sympy.Symbol('sc2')
sp=sympy.Symbol('sp')
oc1=sympy.Symbol('oc1')
oc2=sympy.Symbol('oc2')
op=sympy.Symbol('op')
wc1=sympy.Symbol('wc1')
wc2=sympy.Symbol('wc2')
wp=sympy.Symbol('wp')
E0=sympy.Symbol('E0')
E1=sympy.Symbol('E1')
E2=sympy.Symbol('E2')
E3=sympy.Symbol('E3')
h=sympy.Symbol('h')
t=sympy.Symbol('t')
i = complex(0,1)

H = [[E0, -h/2*oc1*sympy.exp(i*(wc1*t + fic1)), -h/2*op*sympy.exp(i*(wp*t + fip))], [-h/2*oc1*sympy.exp(-i*(wc1*t + fic1)), E1,-h/2*oc2*sympy.exp(i*(wc2*t + fic2))], [-h/2*op*sympy.exp(-i*(wp*t + fip)), -h/2*oc2*sympy.exp(-i*(wc2*t + fic2)), E2]]

Ur1 = [[1, 0, 0], [0, sympy.exp(-i*(wc1*t + fic1)),  0], [0, 0,sympy.exp(-i*(wp*t + fip))]]
Ur2 =  [[1, 0, 0], [0, sympy.exp(i*(wc1*t + fic1)),  0], [0, 0,sympy.exp(i*(wp*t + fip))]]

Urd = [[0, 0, 0], [0, i*(wc1)*sympy.exp(i*(wc1*t + fic1)), 0], [0, 0,i*(wp)*sympy.exp(i*(wp*t + fip))]]

result0 = [[0,0,0],
         [0,0,0],
         [0,0,0]]
result1 = [[0,0,0],
         [0,0,0],
         [0,0,0]]
result2 = [[0,0,0],
         [0,0,0],
         [0,0,0]]

def Ur2H():
    for i in range(len(Ur2)):
       for j in range(len(H[0])):
           for k in range(len(H)):
               result0[i][j] += Ur2[i][k] *H[k][j]
    b = result0
#    print(b)
    return b

U2H = Ur2H()

def Ur2HUr1():
    for i in range(len(U2H)):
        for j in range (len(Ur1[0])):
            for k in range (len(Ur1)):
                result1[i][j]+=U2H[i][k]*Ur1[k][j]
    a = result1
    for r in a:
        print(r)
    return a

U2HU1 = Ur2HUr1()

def UrdUr1():
    for i in range(len(Urd)):
        for j in range (len(Ur1[0])):
            for k in range (len(Ur1)):
                result2[i][j]+=Urd[i][k]*Ur1[k][j]*i*h
    a = result2
    print (a)
    return a 
UdU1= UrdUr1()

rezultat=[[0,0,0],
         [0,0,0],
         [0,0,0]]
for i in range(len(U2HU1)):
   for j in range(len(U2HU1[0])):
       rezultat[i][j] = U2HU1[i][j] + UdU1[i][j]

for r in rezultat:
    print(r)


