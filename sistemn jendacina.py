# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:06 2019

@author: nina
"""


import numpy as np
import scipy.sparse.linalg
import scipy.sparse

dp=0.5
dc=0
op=0.1
oc=1
ii=complex(0,1)
g00=0
g01=0
g02=0
g10=0
g11=0
g12=0
g20=0
g21=0
g22=0
G00=0
G01=0
G02=0

a=[[0,0,-1/2*op*ii+G01,0,0,0,1/2*op*ii,0,G02],
           [0,ii*(dp-dc)-g01,-1/2*oc*ii,0,0,0,0,1/2*op*ii,0],
           [-1/2*op*ii,-1/2*oc*ii,-ii*dp-g02,0,0,0,0,0,1/2*op*ii],
           [0,0,0,-ii*(dp-dc)-g10,-1/2*op*ii,0,1/2*oc*ii,0,0],
           [0,0,0,0,-G01,-1/2*oc*ii,0,1/2*oc*ii,0],
           [0,0,0,-1/2*op*ii,1/2*oc*ii,-ii*((dp-dc)+dp)-g12,0,0,1/2*oc*ii],
           [1/2*op*ii,0,0,1/2*oc*ii,0,0,ii*dp-g20,0,1/2*op*ii],
           [0,1/2*op*ii,0,0,1/2*oc*ii,0,0,ii*((dp-dc)-dp)-g21,1/2*oc*ii],
           [0,0,1/2*op*ii,0,0,1/2*oc*ii,-1/2*op*ii,-1/2*oc*ii,-G02]]

A = scipy.sparse.csr_matrix(a, dtype=np.cfloat)

b=[1,1,1,1,1,1,1,1,1]
print(np.linalg.solve(a,b))
