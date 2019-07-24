# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 22:50:03 2019

@author: nina
"""
import numpy as np
i = complex(0,1)
A = [[3*i,3],[0,i]]
b = [0,0]
print(np.linalg.solve(A,b))
