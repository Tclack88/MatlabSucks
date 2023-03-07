from math import sin, cos, pi
from numpy import array, matmul as mul


R1 = array([[cos(pi/6),     0,      -sin(pi/6)],
            [0,             1,      0],
            [sin(pi/6),    0,      cos(pi/6)]]) # R0-> 1

R2 = array([[1,         0,              0],
            [0,         cos(pi/6),     -sin(pi/6)],
            [0,         sin(pi/6),    cos(pi/6)]]) #R1->2


r1 = array([1,2,3])
r2 = array([3,2,1])
r3 = array([2,1,3])


a1 = mul(R1,r1)
a2 = a1+r2
a3 = mul(R2,a2)
a4 = a3+r3

a5 = mul(R2.T,a4)
a6 = mul(R1.T,a5)

print(a6)

a7 = mul(R1,r1) +
