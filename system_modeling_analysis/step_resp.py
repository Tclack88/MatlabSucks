import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import mpmath as mp
from math import sqrt

x,s = sym.symbols('x,s')
# Stuff that changes (time range, function)
t = np.linspace(0.4,1,100)
#function = 16/(s**2+16*s+16)
K,w,Z = 1,sqrt(50),0
def second_order_sys(K,Z,w,s):
    return K*s**2/(s**2+2*Z*w*s+w**2)

def second_order_sys2(K,Z,w,s):
    return K*w**2/(s**2+2*Z*w*s+w**2)

function = second_order_sys2(K,Z,w,s)

# symbolic solution to get equation if necessary
#print the output of this if it's complicated to see what goes in f(s) below
step_resp = sym.simplify(function*(1/s),rational=True)

print('step_resp')
print(step_resp)
print('inverse laplace:') # doesn't work well with complicated equations
print(sym.inverse_laplace_transform(step_resp,s,x))

#import sys; sys.exit(0)

def f(s): # the step response (printed above for convenience) as a function for mpmath
    return mp.cos(5*sqrt(2)*s)#*Heaviside(s)

#mathematical solution and plot:
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
glist = list(G)
c = np.array(glist, dtype=float)
print(c[0],c[-1]) # get initial and final (steady-state) values
plt.plot(t,c)
plt.show()


