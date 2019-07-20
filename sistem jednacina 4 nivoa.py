# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 00:36:36 2019

@author: nina
"""
import sympy
sympy.init_printing()

i = complex(0,1)
ro00=sympy.Symbol('ro00')
ro01=sympy.Symbol('ro01')
ro02=sympy.Symbol('ro02')
ro10=sympy.Symbol('ro10')
ro11=sympy.Symbol('ro11')
ro12=sympy.Symbol('ro12')
ro20=sympy.Symbol('ro20')
ro21=sympy.Symbol('ro21')
ro22=sympy.Symbol('ro22')
ro30=sympy.Symbol('ro30')
ro31=sympy.Symbol('ro31')
ro32=sympy.Symbol('ro32')
ro33=sympy.Symbol('ro33')
ro03=sympy.Symbol('ro03')
ro13=sympy.Symbol('ro13')
ro23=sympy.Symbol('ro23')

dp=sympy.Symbol('dp')
dc1=0
dc2=0
oc1=1
oc2=1
op=0.1
h=sympy.Symbol('h')

M=sympy.Matrix([[0,-i/2*oc1,0,-i/2*op,i/2*oc1,1,0,0,i/2*oc2,0,1,0,i/2*op,0,0,1],
 [0,0,i/2*oc2,0,0,0,0,0,-i/2*oc2,0,1,0,0,0,0,0],
 [0,i/2*oc1,0,0,-i/2*oc1,1,0,0,0,0,0,0,0,0,0,0],
 [0,0,0,i/2*op,0,0,0,0,0,0,0,0,-i/2*op,0,0,1],
 [-i/2*oc1,1+i*dc1,0,0,0,i/2*oc1,0,0,0,i/2*oc2,0,0,0,i/2*op,0,0],
 [-i/2*oc2,0,i*dc2-1,0,0,0,i/2*oc1,0,0,0,i/2*oc2,0,0,0,i/2*op,0],
 [-i/2*op,0,0,i*dp-1,0,0,0,i/2*oc1,0,0,0,i/2*oc2,0,0,0,i/2*op],
 [0,0,i/2*oc1,0,-i/2*oc2,0,i*(dc2-dc1)-1,0,0,0,0,0,0,0,0,0],
 [0,0,0,i/2*oc1,-i/2*op,0,0,i*(dp-dc1)-1,0,0,0,0,0,0,0,0],
 [0,0,0,i/2*oc2,0,0,0,0,-i/2*op,0,0,i*(dp-dc2)-1,0,0,0,0],
 [i/2*oc1,0,0,0,i*dc1-1,-i/2*oc1,-i/2*oc2,-i/2*op,0,0,0,0,0,0,0,0],
 [i/2*oc2,0,0,0,0,0,0,0,i*dc2-1,-i/2*oc1,i/2*oc2,i/2*op,0,0,0,0],
 [0,i/2*oc2,0,0,0,0,0,0,i/2*oc1,i*(dc1+dc2)-1,0,0,0,0,0,0],
 [i/2*op,0,0,0,0,0,0,0,0,0,0,0,-i*dp-1,-i/2*oc1,-i/2*oc2,-i/2*op],
 [0,0,i/2*op,0,0,0,0,0,0,0,0,0,-i/2*oc2,-1,i*(dp-dc2),0],
 [0,0,i/2*op,0,0,0,0,0,0,0,0,0,-i/2*oc2,0,i*(dp-dc2)-1,0]])


N = sympy.Matrix([ro00,ro01,ro02,ro03,ro10,ro11,ro12,ro13,ro20,ro21,ro22,ro23,ro30,ro31,ro32,ro33])
K = sympy.Matrix([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]) 
#print(len(M[0]),len(M[1]),len(M[2]),len(M[3]),len(M[4]),len(M[5]),len(M[6]),len(M[7]),len(M[8]),len(M[9]),len(M[10]),len(M[11]),len(M[12]),len(M[13]),len(M[14]),len(M[15]))
a = sympy.linsolve((M,N),K)
print(a)
