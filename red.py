# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:26:08 2019

@author: nina
"""
import sympy

dc=0.
dp=0.
oc=0.5
op=0.5

ro00=sympy.Symbol('ro00')
ro01=sympy.Symbol('ro01')
ro02=sympy.Symbol('ro02')
ro10=sympy.Symbol('ro10')
ro11=sympy.Symbol('ro11')
ro12=sympy.Symbol('ro12')
ro20=sympy.Symbol('ro20')
ro21=sympy.Symbol('ro21')
ro22=sympy.Symbol('ro22')
Gamma00=1.
Gamma01=1.
Gamma02=1.
Gamma21=1.
Gamma20=1.
Gamma10=1.
Gamma11=1.
Gamma12=1.
Gamma22=1.
I = complex(0,1)

gamma00=1.
gamma01=1.
gamma02=1.
gamma10=1.
gamma11=1.
gamma12=1.
gamma20=1.
gamma21=1.
gamma22=1.
l = 0.2

ro = [ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22]
#lambda na 0
def na0(ro00,ro01,r02,ro10,ro11,ro12,ro20,ro21,ro22):
    f1 = Gamma01*ro11 + Gamma02*ro22 - I*0.5*op*ro20
    f2 = -gamma01*ro01 - 1.0*I*(0.5*oc*ro02 - ro01*(dc - dp))
    f3 = -gamma02*ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01)
    f4=-gamma10*ro10 - 1.0*I*(-0.5*oc*ro20 + ro10*(dc - dp))
    f5=-Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21)
    f6=-gamma12*ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + ro12*(dc - dp))
    f7=-gamma20*ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10)
    f8=-gamma21*ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - ro21*(dc - dp))
    f9=-Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21)
    return(f1,f2,f3,f4,f5,f6,f7,f8,f9)
    
#lambda na 1
def na1(ro00,ro01,r02,ro10,ro11,ro12,ro20,ro21,ro22):
    f1 = (Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*op*ro02 - 0.5*op*ro20))*l
    f2 = (-gamma01*ro01 - 1.0*I*(0.5*oc*ro02 - 0.5*op*ro21 - ro01*(dc - dp)))*l
    f3 = (-gamma02*ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01 + 0.5*op*ro00 - 0.5*op*ro22))*l
    f4 = (-gamma10*ro10 - 1.0*I*(-0.5*oc*ro20 + 0.5*op*ro12 + ro10*(dc - dp)))*l
    f5 = (-Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21))*l
    f6 = (-gamma12*ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + 0.5*op*ro10 + ro12*(dc - dp)))*l
    f7 = (-gamma20*ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10 - 0.5*op*ro00 + 0.5*op*ro22))*l
    f8 = (-gamma21*ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - 0.5*op*ro01 - ro21*(dc - dp)))*l
    f9 = (-Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21 - 0.5*op*ro02 + 0.5*op*ro20))*l
    return(f1,f2,f3,f4,f5,f6,f7,f8,f9)

#lambda na 2
def na2(ro00,ro01,r02,ro10,ro11,ro12,ro20,ro21,ro22):
    f1 = (Gamma01*ro11 + Gamma02*ro22 - 1.0*I*(0.5*op*ro02 - 0.5*op*ro20))*l**2
    f2 = (-gamma01*ro01 - 1.0*I*(0.5*oc*ro02 - 0.5*op*ro21 - ro01*(dc - dp)))*l**2
    f3 = (-gamma02*ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01 + 0.5*op*ro00 - 0.5*op*ro22))*l**2
    f4 = (-gamma10*ro10 - 1.0*I*(-0.5*oc*ro20 + 0.5*op*ro12 + ro10*(dc - dp)))*l**2
    f5 = (-Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21))*l**2
    f6 = (-gamma12*ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + 0.5*op*ro10 + ro12*(dc - dp)))*l**2
    f7 = (-gamma20*ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10 - 0.5*op*ro00 + 0.5*op*ro22))*l**2
    f8 = (-gamma21*ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - 0.5*op*ro01 - ro21*(dc - dp)))*l**2
    f9 = (-Gamma02*ro22 - Gamma12*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21 - 0.5*op*ro02 + 0.5*op*ro20))*l**2
    return(f1,f2,f3,f4,f5,f6,f7,f8,f9)
#lambda na 3
def na3(ro00,ro01,r02,ro10,ro11,ro12,ro20,ro21,ro22):
    g1 = (-1.0*I*(0.5*op*ro02 - 0.5*op*ro20))*l**3
    g2 = (1.0*I*(0.5*op*ro21))*l**3
    g3 = (1.0*I*(-0.5*op*ro00 + 0.5*op*ro22))*l**3
    g4 = (-0.5*I*op*ro12)*l**3
    #g5 = (-Gamma01*ro11 + Gamma12*ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21))*l
    g6 = (- 1.0*I*(0.5*op*ro10))*l**3
    g7 = (- 1.0*I*(- 0.5*op*ro00 + 0.5*op*ro22))*l**3
    g8 = (- 1.0*I*(- 0.5*op*ro01))*l**3
    g9 = (- 1.0*I*(- 0.5*op*ro02 + 0.5*op*ro20))*l**3
    return(g1,g2,g3,g4,g6,g7,g8,g9)

nultiStepen = na0(ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22)
prviStepen = na1(ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22)
drugiStepen = na2(ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22)
#treciStepen = na3(ro)

jednacina1 = sympy.solve(prviStepen,[ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22])
jednacina2 = sympy.solve(drugiStepen,[ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22])
jednacina0 = sympy.solve(nultiStepen,[ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22])
print(jednacina0)