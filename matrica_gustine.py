# -*- coding: utf-8 -*-
"""
Created on Sat Apr  6 01:11:06 2019

@author: nina
"""
#

import sympy
import numpy
from sympy import *
from sympy.solvers.solvers import solve_linear_system
import scipy.sparse.linalg
import scipy.sparse
import matplotlib.pyplot as plt
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
ro10=sympy.Symbol('ro10')
ro11=sympy.Symbol('ro11')
ro12=sympy.Symbol('ro12')
ro20=sympy.Symbol('ro20')
ro21=sympy.Symbol('ro21')
ro22=sympy.Symbol('ro22')
Gamma00=sympy.Symbol('Gamma00')
Gamma01=sympy.Symbol('Gamma01')
Gamma02=sympy.Symbol('Gamma02')
Gamma21=sympy.Symbol('Gamma21')
Gamma20=sympy.Symbol('Gamma20')
Gamma10=sympy.Symbol('Gamma10')
Gamma11=sympy.Symbol('Gamma11')
Gamma12=sympy.Symbol('Gamma12')
Gamma22=sympy.Symbol('Gamma22')


gamma00=sympy.Symbol('gamma00')
gamma01=sympy.Symbol('gamma01')
gamma02=sympy.Symbol('gamma02')
gamma10=sympy.Symbol('gamma10')
gamma11=sympy.Symbol('gamma11')
gamma12=sympy.Symbol('gamma12')
gamma20=sympy.Symbol('gamma20')
gamma21=sympy.Symbol('gamma21')
gamma22=sympy.Symbol('gamma22')
h = sympy.Symbol('h')
ro=[[ro00, ro01, ro02],[ro10, ro11, ro12],[ro20, ro21, ro22]]
gamma=[[gamma00, gamma01, gamma02],[gamma10, gamma11, gamma12],[gamma20, gamma21, gamma22]]
Gamma=[[Gamma00, Gamma01, Gamma02],[Gamma10, Gamma11, Gamma12],[Gamma20, Gamma21, Gamma22]]
ii = sympy.I

result0 = [[0,0,0],
         [0,0,0],
         [0,0,0]]
result1 = [[0,0,0],
         [0,0,0],
         [0,0,0]]
result = [[0,0,0],
         [0,0,0],
         [0,0,0]]

brojac=[]
imaginarnaGustina21=[]
realnaGustina21=[]
imaginarnaGustina02=[]
realnaGustina02=[]
imaginarnaGustina10=[]
realnaGustina10=[]
H=[[0, 0, -1/2*op],[0, -(dp-dc), -1/2*oc],[-1/2*op, -1/2*oc, -dp]]
for i in range(len(H)):
   for j in range(len(ro[0])):
       for k in range(len(ro)):
           result0[i][j] += H[i][k] * ro[k][j]

for i in range(len(ro)):
   for j in range(len(H[0])):
       for k in range(len(H)):
           result1[i][j] += ro[i][k] * H[k][j]
           
for i in range (3):
    for j in range (3):
        result[i][j]=-complex(0,1)*(result0[i][j]-result1[i][j])
        if i!=j:
            result[i][j]+=-gamma[i][j]*ro[i][j]
        else:
            if i == 0:
                result[i][j]=result[i][j]+Gamma[0][1]*ro[1][1]+Gamma[0][2]*ro[2][2]
            elif i == 1:
                result[i][j]=result[i][j]-Gamma[0][1]*ro[1][1]+Gamma[1][2]*ro[2][2]
            else:
                result[i][j]=result[i][j]-Gamma[0][2]*ro[2][2]-Gamma[1][2]*ro[2][2]
#print(result)                
for i in range (3):
    for j in range (3):
        locals()["jna"+str(i)+str(j)] = result[i][j]
        print(result[i][j])                

jna = ro00 + ro11 + ro22 - 1
#jednacina = sympy.solve([jna01,jna02,jna10,jna11,jna12,jna20,jna21,jna00,jna],[ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22])

##ovo je ako zelim da nadjem konkretnu vrednost za ro odn da resim jnu
#for i in range (100):
#    dp = (i-50)/50
#    print(dp)
#    H=[[0, 0, -1/2*op*h],[0, -(dp-dc)*h, -1/2*op*h],[-1/2*op*h, -1/2*oc*h, -dp*h]]
#    for i in range(len(H)):
#       for j in range(len(ro[0])):
#           for k in range(len(ro)):
#               result0[i][j] += H[i][k] * ro[k][j]
#    
#    for i in range(len(ro)):
#       for j in range(len(H[0])):
#           for k in range(len(H)):
#               result1[i][j] += ro[i][k] * H[k][j]
#               
#    for i in range (3):
#        for j in range (3):
#            result[i][j]=-ii*(result0[i][j]-result1[i][j])/h
#            if i!=j:
#                result[i][j]+=-gamma[i][j]*ro[i][j]
#            else:
#                #jedina stvar za koju nismo sigurne u kodu je da li smo dobro protumacile
#                #kako se odredjuju jednacine na dijagonali 
#                if i == 0:
#                    result[i][j]=result[i][j]+Gamma[0][1]*ro[1][1]+Gamma[0][2]*ro[2][2]
#                elif i == 1:
#                    result[i][j]=result[i][j]-Gamma[0][1]*ro[1][1]+Gamma[1][2]*ro[2][2]
#                else:
#                    result[i][j]=result[i][j]-Gamma[0][2]*ro[2][2]-Gamma[1][2]*ro[2][2]
#    print(result)                
#    for i in range (3):
#        for j in range (3):
#            locals()["jna"+str(i)+str(j)] = result[i][j]
##            print(result[i][j])
#    jna = ro00 + ro11 + ro22 - 1
#    jednacina = sympy.solve([jna01,jna02,jna10,jna11,jna12,jna20,jna21,jna00,jna],[ro00,ro01,ro02,ro10,ro11,ro12,ro20,ro21,ro22])
#    
#    brojac.append(dp)
#    imaginarnaGustina21.append(im(jednacina[ro21]))
#    realnaGustina21.append(re(jednacina[ro21]))
#    imaginarnaGustina02.append(im(jednacina[ro02]))
#    realnaGustina02.append(re(jednacina[ro02]))
#    imaginarnaGustina10.append(im(jednacina[ro10]))
#    realnaGustina10.append(re(jednacina[ro10]))
##print(realnaGustina)
##print (brojac)
#plt.plot(brojac,realnaGustina21)
#plt.title("Realni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho21")
#plt.show()
#plt.title("Imaginarni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho21")
#plt.plot(brojac,imaginarnaGustina21)
#plt.show()    
#plt.plot(brojac,realnaGustina10)
#plt.title("Realni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho10")
#plt.show()
#plt.title("Imaginarni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho10")
#plt.plot(brojac,imaginarnaGustina10)
#plt.show()   
#plt.plot(brojac,realnaGustina02)
#plt.title("Realni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho02")
#plt.show()
#plt.title("Imaginarni deo")
#plt.xlabel("delta p")
#plt.ylabel("rho02")
#plt.plot(brojac,imaginarnaGustina02)
#plt.show()   

#ovo je ako zwwlim da nadjem resenje u opstim brojevima
