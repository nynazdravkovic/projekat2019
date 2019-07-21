# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 10:23:48 2019

@author: nina
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Jul 21 01:08:44 2019

@author: nina
"""

from IPython.display import display
import sympy
from sympy import Function, sqrt, dsolve, Eq, Derivative
from sympy import solve, Poly, Eq, Function, exp
from sympy import Indexed, IndexedBase, Tuple, sqrt


t = sympy.Symbol('t')
d=0
o=0.5
ro00=Function('ro00')(t)
ro01=Function('ro01')(t)
ro10=Function('ro10')(t)
ro11=Function('ro11')(t)
I = sympy.Symbol('I')

eq1 = Eq(Derivative(ro00,t), ro11 - I*(o*ro01/2 - o*ro10/2))
eq2 = Eq(Derivative(ro01,t), -ro01 - I*(-d*ro01 + o*ro00/2 - o*ro11/2))
eq3 = Eq(Derivative(ro10,t),-ro10 - I*(d*ro10 - o*ro00/2 + o*ro11/2))
eq4 = Eq(Derivative(ro11,t), -ro11 - I*(-o*ro01/2 + o*ro10/2))

"resenje je sa konstantama C1, C2 i C3"
soln = dsolve((eq1, eq2, eq3, eq4),ics={ro00(0):1, ro01(0):0, ro10(0):0, ro11(0):0})
display(soln)
display(soln)
print(soln)