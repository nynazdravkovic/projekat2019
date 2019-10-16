

import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
import scipy.integrate as integrate
from scipy.fftpack import fft, ifft
#W konfiguracija
fi = 0
oc0 = 1.
oc1 = 1.
oc2 = 1.
d0 = 0.
d1 = 0.
d2 = 0.
dm = 0.
op = 0.5
W = 1.
e = math.exp
Gamma01 = 1.
Gamma02 = 1.
Gamma12 = 1.
Gamma03 = 1.
Gamma13 = 1.
Gamma23 = 1.
gamma23 = 1.
gamma0 = 1.
gamma02 = 1.
gamma1 = 1.
gamma12 = 1.
gamma2 = 1.
gamma21 = 1.
gamma03 = 1.
gamma3 = 1.
gamma13 = 1.
gamma30 = 1.
gamma31 = 1.
gamma32 = 1.
gamma33 = 1.
fi = 0
oc0 = 1.
d0 = 0.
d1 = 0.
dm = 0.
op = 0.5
I = complex(0,1)
n=10
c=1
udaljenost = np.linspace(0,5,n)
vreme = np.linspace(0,5,n)
E = np.empty((n,n))
epsilon = 1

#vraca ro02 sa tildom
def f (w):
    nultiStepen=[[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], #00
                [-I*oc0,I*d0-gamma0,-I*oc1,0,0,I*oc0,0,0,0,0,0,0,0,0,0,0], #01
                 [0,-I*oc1,I*(d0+d1)-gamma0-gamma1,-I*oc2,0,0,-I*oc0,0,0,0,0,0,0,0,0,0], #02
                 [0, 0, -I*oc2, I*(d0+dm)-gamma0-gamma2-gamma3,0,0,0,I*oc0,0,0,0,0,0,0,0,0], #03
                 [I*oc0, 0,0,0,I*d0-gamma0,-I*oc0,0,0, I*oc0,0,0,0,0,0,0,0], #10
                 [0,-I*oc0,0,0,I*oc0,-I*gamma0,-I*oc1,0,0,I*oc1,-gamma1,0,0,0,0,-gamma3],#11
                 [0,0,I*oc0,0,0,-I*oc1, -I*(d0+d1)-gamma1,-I*oc2,0,0,I*oc1,0,0,0,0,0], #12
                 [0,0,0,I*oc0,0,0,-I*oc2,-I*dm-gamma2-gamma3,0,0,0,I*oc1,0,0,0,0], #13
                 [0,0,0,0,I*oc1,0,0,0, -I*(d0+d1)-gamma0-gamma1,-I*oc0,0,0,I*oc2,0,0,0], #20
                 [0,0,0,0,0,I*oc2,0,0,-I*oc0, I*(d0+d1)-gamma1, -I*oc1,0,0,I*oc2,0,0], #21
                 [0,0,0,0,0,0,I*oc1,0,0,-I*oc1,-gamma1,-I*oc2,0,0,I*oc2,-gamma2], #22
                 [0,0,0,0,0,0,0,I*oc1,0,0,-I*oc2, -I*(d0+d1+d2)-gamma2-gamma3,0,I*oc1,0,-I*oc2], #23
                 [0,0,0,0,0,0,0,0,I*oc2,0,0,0, -I*(d0+dm)-gamma0-gamma2-gamma3,-I*oc0,0,0], #30
                 [0,0,0,0,0,0,0 ,0 ,0, I*oc2,0,0,- I*oc0, I*dm-gamma2-gamma3, -I*oc1,0],
                 [0,0,0,0,0,0,0,0,0,0,I*oc2,0,0,-I*oc1, I*(d0+d1+d2)-gamma2-gamma3,-I*oc2],
                 [0,0,0,0,0,0,0,0,0,0,0,I*oc2,0,0,-I*oc2,-gamma3-gamma2]
                 ]
    
    
    nule = np.zeros(16,dtype='complex')
    nule[0] = 1
    resenje1 = np.linalg.solve(nultiStepen,nule)
# ovde resavam jednacine za prvi stepen a ono sto stoji uz op uzimam iz 
#prethodne matrice jer uz njega stoji nesto sto mnoz lambda na 0 :)
    prviStepen=[[1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1], #00
                [-I*oc0,I*d0-gamma0,-I*oc1,0,0,I*oc0,0,0,0,0,0,0,0,0,0,0], #01
                 [0,-I*oc1,I*(d0+d1)-gamma0-gamma1,-I*oc2,0,0,-I*oc0,0,0,0,0,0,0,0,0,0], #02
                 [0, 0, -I*oc2, I*(d0+dm)-gamma0-gamma2-gamma3,0,0,0,I*oc0,0,0,0,0,0,0,0,0], #03
                 [I*oc0, 0,0,0,I*d0-gamma0,-I*oc0,0,0, I*oc0,0,0,0,0,0,0,0], #10
                 [0,-I*oc0,0,0,I*oc0,-I*gamma0,-I*oc1,0,0,I*oc1,-gamma1,0,0,0,0,-gamma3],#11
                 [0,0,I*oc0,0,0,-I*oc1, -I*(d0+d1)-gamma1,-I*oc2,0,0,I*oc1,0,0,0,0,0], #12
                 [0,0,0,I*oc0,0,0,-I*oc2,-I*dm-gamma2-gamma3,0,0,0,I*oc1,0,0,0,0], #13
                 [0,0,0,0,I*oc1,0,0,0, -I*(d0+d1)-gamma0-gamma1,-I*oc0,0,0,I*oc2,0,0,0], #20
                 [0,0,0,0,0,I*oc2,0,0,-I*oc0, I*(d0+d1)-gamma1, -I*oc1,0,0,I*oc2,0,0], #21
                 [0,0,0,0,0,0,I*oc1,0,0,-I*oc1,-gamma1,-I*oc2,0,0,I*oc2,-gamma2], #22
                 [0,0,0,0,0,0,0,I*oc1,0,0,-I*oc2, -I*(d0+d1+d2)-gamma2-gamma3,0,I*oc1,0,-I*oc2], #23
                 [0,0,0,0,0,0,0,0,I*oc2,0,0,0, -I*(d0+dm)-gamma0-gamma2-gamma3,-I*oc0,0,0], #30
                 [0,0,0,0,0,0,0 ,0 ,0, I*oc2,0,0,- I*oc0, I*dm-gamma2-gamma3, -I*oc1,0],
                 [0,0,0,0,0,0,0,0,0,0,I*oc2,0,0,-I*oc1, I*(d0+d1+d2)-gamma2-gamma3,-I*oc2],
                 [0,0,0,0,0,0,0,0,0,0,0,I*oc2,0,0,-I*oc2,-gamma3-gamma2]
                 ]
    resenjePrvog = [1,
                    I*op*resenje1[3],
                    0,
                    I*op*(resenje1[1]),
                    -I*op*resenje1[12],
                    I*op*(resenje1[7]-resenje1[13]),
                    -I*op*resenje1[14],
                    I*op*(resenje1[5]-resenje1[15]),
                    0,
                    I*op*(resenje1[7]),
                    0,
                    I*op*resenje1[9],
                    -I*op*resenje1[4],
                    I*op*(resenje1[15]-resenje1[5]),
                    -I*op*resenje1[6],
                    I*op*(resenje1[13]-resenje1[7])]
    resenje2 = np.linalg.solve(prviStepen,resenjePrvog)
    return(resenje2[13]/(I*op))

#funkcija koja racuna povrsinu ispod funkcije za odvojei imaginarni i realni deo
def integral(func,a,b):
    def realnaFja(x):
        return scipy.real(func(x))
    def imaginarnaFja(x):
        return scipy.imag(func(x))
    integralRe = integrate.quad(realnaFja,a,b)
    integralIm = integrate.quad(imaginarnaFja,a,b)
    return (integralRe[0]+I*integralIm[0])

#ovo smo izracunale rucno 
def gaus(w,W):
    a = epsilon*W/np.sqrt(2*np.log(2))*np.exp(-(w)**2*(W)**2/(8*np.log(8)))/np.sqrt(2)
    return(a)

#Ep0 je niz koji se odmah racuna, to je gaus koji je prethodni prosao kroz fft 
#da bi bio gaus treba da t ide u neg vrednosti 
#k je mesto u nizu a t je vreme zato sto je vreme i negaivno a niz je izracunat tim vremenop pre
def ispodIntegrala(w,t,k, z):
    return(gaus(w,W)*np.exp(-I*w*(t-z/c)-f(w)*z))
    
def resavanje(distance,time,W):
    matricaResenjaRe = np.zeros((len(time),len(distance)))
    matricaResenjaIm = np.zeros((len(time),len(distance)))
    matricaModuo=np.zeros((len(time),len(distance)))
    for  i in range(len(time)):
        t = time[i]
        for j in range(len(distance)):
            z = distance[j]
            Ep = integral(lambda w: ispodIntegrala(w,t,i,z),-np.inf,np.inf)
#            print(Ep)
            matricaResenjaRe[i][j] = 1/(np.sqrt(2*np.pi))*Ep.real
            matricaResenjaIm[i][j] = 1/(np.sqrt(2*np.pi))*Ep.imag
            matricaModuo[i][j]=np.sqrt(matricaResenjaRe[i][j]**2+matricaResenjaIm[i][j]**2)
    return(matricaResenjaRe,matricaResenjaIm,matricaModuo)
    
polje = resavanje(udaljenost,vreme,W)
poljeIm = polje[1]
poljeRe = polje[0]
poljeModuo = polje[2]

#np.savetxt('Ep1.csv',polje, delimiter=',')
N=n*n
T = np.meshgrid(vreme,udaljenost)
#plt.contourf(T[1],T[0],poljeRe)
plt.contourf(T[1], T[0], poljeRe, N, cmap="viridis")
plt.colorbar()

plt.show()
plt.contourf(T[1], T[0], poljeIm, N, cmap="viridis")
plt.colorbar()
plt.show()
plt.contourf(T[1], T[0], poljeModuo, N, cmap="viridis")
plt.colorbar()
plt.title("Moduo probnog polja, Ep")
plt.xlabel("x[m]")
plt.ylabel("t[s]")
plt.show()
#
#
