# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 17:37:17 2019

@author: nina
"""
"kod za resavanje sistema ODE sa razdvajanjem promenljivih"
from sympy import init_printing
from IPython.display import display
from sympy import *
from sympy.abc import x, y, z, w
import sympy
#integrate(t)
#u ovom kodu resavam dif jne 
from sympy import Function, Indexed, Tuple, sqrt, dsolve, solve, Eq, Derivative, sin, cos, symbols
from sympy import solve, Poly, Eq, Function, exp
"zbog biblioteke sympy.abc rho00=x rho01=y rho10=z i rho11=w"

f = Function('f')
dc=sympy.Symbol('dc')
dp=sympy.Symbol('dp')
oc=sympy.Symbol('oc')
op=sympy.Symbol('op')
I=sympy.Symbol('I')

ro00=sympy.Symbol('ro00', cls = Function, Function = True)
ro01=sympy.Symbol('ro01', cls = Function, Function = True)
ro02=sympy.Symbol('ro02', cls = Function, Function = True)
ro10=sympy.Symbol('ro10', cls = Function, Function = True)
ro11=sympy.Symbol('ro11', cls = Function, Function = True)
ro12=sympy.Symbol('ro12', cls = Function, Function = True)
ro20=sympy.Symbol('ro20', cls = Function, Function = True)
ro21=sympy.Symbol('ro21', cls = Function, Function = True)
ro22=sympy.Symbol('ro22', cls = Function, Function = True)
#x, y, z, w = symbols("x y z w", cls = Function, Function = True)
#eq1 = Eq(Derivative(y(t),t), i*(-d*h*y(t) + h*o*x(t)/2 - h*o*w(t)/2)/h)
#eq2 = Eq(Derivative(z(t),t), i*(d*h*z(t) - h*o*x(t)/2 + h*o*w(t)/2)/h)
#eq3 = Eq(Derivative(w(t),t),i*(-h*o*y(t)/2 + h*o*z(t)/2)/h )
eq1 = Eq(Derivative(ro00(t),t),ro11 + ro22 - 1.0*I*(0.5*op*ro02 - 0.5*op*ro20))
eq2 = Eq(Derivative(ro01(t),t),-ro01 - 1.0*I*(0.5*oc*ro02 - 0.5*op*ro21 - ro01*(dc - dp)))
eq3 = Eq(Derivative(ro02(t),t),-ro02 - 1.0*I*(dp*ro02 + 0.5*oc*ro01 + 0.5*op*ro00 - 0.5*op*ro22))
eq4 = Eq(Derivative(ro10(t),t),-ro10 - 1.0*I*(-0.5*oc*ro20 + 0.5*op*ro12 + ro10*(dc - dp)))
eq5 = Eq(Derivative(ro11(t),t),-ro11 + ro22 - 1.0*I*(0.5*oc*ro12 - 0.5*oc*ro21))
eq6 = Eq(Derivative(ro12(t),t),-ro12 - 1.0*I*(dp*ro12 + 0.5*oc*ro11 - 0.5*oc*ro22 + 0.5*op*ro10 + ro12*(dc - dp)))
eq7 = Eq(Derivative(ro20(t),t),-ro20 - 1.0*I*(-dp*ro20 - 0.5*oc*ro10 - 0.5*op*ro00 + 0.5*op*ro22))
eq8 = Eq(Derivative(ro21(t),t),-ro21 - 1.0*I*(-dp*ro21 - 0.5*oc*ro11 + 0.5*oc*ro22 - 0.5*op*ro01 - ro21*(dc - dp)))
eq9 = Eq(Derivative(ro22(t),t),-2*ro22 - 1.0*I*(-0.5*oc*ro12 + 0.5*oc*ro21 - 0.5*op*ro02 + 0.5*op*ro20))

# Solutions with undetermined constants of integration, C1, C2 and C3:
soln = dsolve((eq1, eq2, eq3, eq4, eq5, eq6, eq7, eq8, eq9))
display(soln)
print(soln)

