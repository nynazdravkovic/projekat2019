# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 21:26:08 2019

@author: nina
"""
import numpy as np

dc=0.
dp=0.
oc=0.5
op=0.5

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

#ro = [ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22]
#lambda na 0   
sistem1 = np.array([[0,0,0,0,Gamma01,0,0,0,Gamma02],
           [0,-gamma01+I*(dc-dp),-I*0.5*oc,0,0,0,0,0,0],
           [0,-I*0.5*oc,-gamma02-I*dp,0,0,0,0,0,0],
           [0,0,0,-gamma10-I*(dc-dp),0,0,I*0.5*oc,0,0],
           [0,0,0,0,-Gamma01,-I*0.5*oc,0,I*0.5*oc,Gamma12],
           [0, 0, 0, 0, -1.0*I*0.5*oc, - gamma12-1.0*I*dc, 0, 0, -1.0*I*(-0.5)*oc], 
           [0, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma20 - 1.0*I*(-dp), 0, 0],
           [0, 0, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma21- 1.0*I*(-dp - (dc - dp)), - 1.0*I*(0.5*oc)],
           [0,0,0,0,0,I*0.5*oc,0,-I*0.5*oc,-Gamma02-Gamma12]])

sistem2 = np.array([[0,0,0,0,Gamma01,0,0,-I*0.5*op,Gamma02],
           [0,-gamma01+I*(dc-dp),-I*0.5*oc,0,0,0,0,-I*0.5*op,0],
           [I*0.5*op,-I*0.5*oc,-gamma02-I*dp,0,0,0,0,0,-I*0.5*op],
           [0,0,0,-gamma10-I*(dc-dp),0,-I*0.5*op,I*0.5*oc,0,0],
           [0,0,0,0,-Gamma01,-I*0.5*oc,0,I*0.5*oc,Gamma12],
           [0, 0, 0, 0, -1.0*I*0.5*oc, - gamma12-1.0*I*(dc-0.5*op), 0, 0, -1.0*I*(-0.5)*oc], 
           [I*0.5*op, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma20 - 1.0*I*(-dp), 0, -I*0.5*op],
           [0, I*0.5*op, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma21- 1.0*I*(-dp - (dc - dp)), - 1.0*I*(0.5*oc)],
           [0,0,I*0.5*op,0,0,I*0.5*oc,-I*0.5*op,-I*0.5*oc,-Gamma02-Gamma12]])

sistem3 = np.array([[0,0,-I*0.5*op,0,Gamma01,0,-*0.5*op,0,Gamma02],
           [0,-gamma01+I*(dc-dp),-I*0.5*oc,0,0,0,0,-I*0.5*op,0],
           [I*0.5*op,-I*0.5*oc,-gamma02-I*dp,0,0,0,0,0,-I*0.5*op],
           [0,0,0,-gamma10-I*(dc-dp),0,-I*0.5*op,I*0.5*oc,0,0],
           [0,0,0,0,-Gamma01,-I*0.5*oc,0,I*0.5*oc,Gamma12],
           [0, 0, 0, 0, -1.0*I*0.5*oc, - gamma12-1.0*I*(dc-0.5*op), 0, 0, -1.0*I*(-0.5)*oc], 
           [I*0.5*op, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma20 - 1.0*I*(-dp), 0, -I*0.5*op],
           [0, I*0.5*op, 0, 0, - 1.0*I*(- 0.5*oc), 0, 0, -gamma21- 1.0*I*(-dp - (dc - dp)), - 1.0*I*(0.5*oc)],
           [0,0,I*0.5*op,0,0,I*0.5*oc,-I*0.5*op,-I*0.5*oc,-Gamma02-Gamma12]])

sistem4 = np.array([[0,0,-I*0.5*op,0,0,0,I*0.5*op,0,0],
           [0,0,0,0,0,0,0,-I*0.5*op,0],
           [I*0.5*op,0,0,0,0,0,0,0,-I*0.5*op],
           [0,0,0,0,0,-I*0.5*op,0,0,0],
           [0,0,0,0,0,0,0,0,0],
           [0,0,0,0,0,I*0.5*op,0,0,0],
           [I*0.5*op,0,0,0,0,0,0,0,-I*0.5*op],
           [0,I*0.5*op,0,0,0,0,0,0,0],
           [0,0,I*0.5*op,0,0,0,-I*0.5*op,0,0]])
levo = np.zeros(9)
print(np.shape(sistem1))
jednacina0 = np.linalg.solve(sistem2,levo)
print(jednacina0)