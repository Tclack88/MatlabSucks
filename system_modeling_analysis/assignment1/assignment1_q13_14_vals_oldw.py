from math import sin, cos, atan, pi, sqrt
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import numpy as np
import sympy as sp

I1,I2 = 1,2
m1,m2=1,2
R,r = 5,3
#D1,D2 = 2,3
g=9.8

c1 = I1+I2 +(m1+m2)*R**2 + m2*(R-r)**2
c2 = 2*m2*R*(R-r)
c3 = I2*(1-R/r) + m2*(R-r)**2
c4 = m2*R*(R-r)

c5 = m2*R*(R-r)
c6 = R*g*(m1+m2)
c7 = m2*g*(R-r)

d1 = I2*(1-R/r) + m2*(R-r)**2
d2 = m2*R*(R-r)
d3 = I2*(1-R/r)**2 + m2*(R-r)**2

d4 = m2*R*(R-r)
d5 = m2*g*(R-r)

#print(f"({c1}+{c2}cos(q2))q1'' + ({c3} + {c4}cos(q2))q2'' - {c5}q2'(2q1'+q2')sin(q2) + {c6}sin(q1) + {c7}sin(q1+q2) + {2}q1' = {2*R}F(t)cos(q1)")
 
#print('\n\n')
#print(f"({d1}+{d2}cos(q2))q1'' + ({d3})q2'' + {d4}q1'^2 * sin(q2) + {d5}sin(q1+q2) + {3}q2' = 0")
print('q13,14')
k1 = 10/((c1+c2*cos(pi))-(c3+c4*cos(pi))*(d1+d2*cos(pi))/d3)

k2 = 10/((c3+c4*cos(pi))-d3*(c1+c2*cos(pi))/(d1+d2*cos(pi)))

print('g3:',k1)
print('g4:',k2)

########### equilibrium #########
print('\n\nq15')
F = 3/2*(2+sqrt(3))*g

def equations(vars):
    q1,q2 = vars
    eq1 = 10*F*cos(q1)+d5*sin(q1+q2)*(c3+c4*cos(q2))/d3 + (c6*sin(q1)+c7*sin(q1+q2))
    eq2 = 10*F*cos(q1)+d5*sin(q1+q2)*(c1+c2*cos(q2))/(d1+d2*cos(q2)) + (c6*sin(q1)+c7*sin(q1+q2))
    return [eq1,eq2]

X,Y = fsolve(equations,(1,1))
X1,Y1 = fsolve(equations,(-1,-1))

print('x:',X)
print('y:', Y)
print()
print(X*180/pi)
print(Y*180/pi)
arg = 2*F/(g*(m1+m2))
theta = 180/pi*atan(arg)
print(theta-180, theta)

############### Jacobian, q18 ########################
# def derivs(vars):
#     q1,q2 = vars
#     D1=-d5/d3*cos(q1+q2)*(c3+c4*cos(q2))-(c6*cos(q1)+c7*cos(q1+q2))+10*F*sin(q1)
#     D2=-d5/d3*(cos(q1+q2)*(c3+c4*cos(q2)) - sin(q1+q2)*(c4*sin(q2)))-c7*cos(q1+q2)
# 
#     D3=(c1+c2*cos(q2))/(d1+d2*cos(q2))*(-d5*cos(q1+q2))-(c6*cos(q1)+c7*cos(q1+q2))+10*F*sin(q1)
#     D4=(d1+d2*cos(q2))*(c2*d5*sin(q2)*cos(q1+q2)+(c1+c2*cos(q2))*(d5*sin(q1+q2)*-d2*sin(q2)))/(d1+d2*cos(q2))**2 - c7*cos(q1+q2)
#     return d1,d2,d3,d4
# 
# q1,q2 = 75*pi/180, 105*pi/180
# D = derivs((q1,q2) )
# A = np.array([[0,0,1,0],
#         [0,     0,      0,  1],
#         [D[0],  D[1],   0,  0],
#         [D[2],  D[3],   0,  0]])
# 
# print(A)

#q1  = q1
#q2  = q2
#q1' = q3
#q2' = q4

print('\n\nQ18')
def check_eqs():
    q2 = 75*pi/180
    q1 = 105*pi/180
    q3 = 0
    q4 = 0
    shit1 = c5*q4*(2*q3+q4)*sin(q2)+ c6*sin(q1) + c7*sin(q1+q2) + 2*q3
    shit2 = d4*q3**2*sin(q2) + d5*sin(q1+q2) + 3*q4
    denom1 = c1+c2*cos(q2) - (c3+c4*cos(q2))*(d1+d2*cos(q2))/d3
    denom2 = c3+c4*cos(q2) - d3*(c1+c2*cos(q2))/(d1+d2*cos(q2))
    eq1 = (10*F*cos(q1)+ shit1 + shit2*(c3+c4*cos(q2))/d3)/denom1
    #eq2 = (10*F*cos(q1) + shit1 + shit2/(d1+d2*cos(q2)))/denom2
    eq2 = (10*F*cos(q1) + shit1 + shit2*(c1+c2*cos(q2))/(d1+d2*cos(q2)))/denom2
    print(eq1)
    print(eq2)

print('checking equations, these should be 0')
check_eqs()
    
    

# Te be continued when above equations fixed

q1,q2,q3,q4 = sp.symbols('q1 q2 q3 q4')
def differentiate(q1,q2,q3,q4):
    #q1 = 75*pi/180
    #q2 = 105*pi/180
    #q3 = 0
    #q4 = 0
    shit1 = c5*q4*(2*q3+q4)*sp.sin(q2)+ c6*sp.sin(q1) + c7*sp.sin(q1+q2) + 2*q3
    shit2 = (d4*q3**2*sp.sin(q2) + d5*sp.sin(q1+q2) + 3*q4)
    denom1 = c1+c2*sp.cos(q2) - (c3+c4*sp.cos(q2))*(d1+d2*sp.cos(q2))/d3
    denom2 = c3+c4*sp.cos(q2) - d3*(c1+c2*sp.cos(q2))/(d1+d2*sp.cos(q2))
    def eq1(q1,q2,q3,q4):
        eq1 = (10*F*sp.cos(q1)+ shit1 + shit2*(c3+c4*sp.cos(q2))/d3)/denom1
        return eq1
    def eq2(q1,q2,q3,q4):
        eq2 = (10*F*sp.cos(q1) + shit1 + shit2*(c1+c2*sp.cos(q2))/(d1+d2*sp.cos(q2)))/denom2
        return eq2
    D1 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq1(q1,q2,q3,q4),q1), "numpy")
    D2 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq1(q1,q2,q3,q4),q2), "numpy")
    D3 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq1(q1,q2,q3,q4),q3), "numpy")
    D4 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq1(q1,q2,q3,q4),q4), "numpy")
    #D2 = sp.diff(eq1(q1,q2,q3,q4),q2)
    #D3 = sp.diff(eq1(q1,q2,q3,q4),q3)
    #D4 = sp.diff(eq1(q1,q2,q3,q4),q4)

    #D5 = sp.diff(eq1(q1,q2,q3,q4),q1)
    #D6 = sp.diff(eq1(q1,q2,q3,q4),q2)
    #D7 = sp.diff(eq1(q1,q2,q3,q4),q3)
    #D8 = sp.diff(eq1(q1,q2,q3,q4),q4)
    D5 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq2(q1,q2,q3,q4),q1), "numpy")
    D6 = sp.lambdify((q1,q2,q3,q4), sp.diff(eq2(q1,q2,q3,q4),q2), "numpy")
    D7= sp.lambdify((q1,q2,q3,q4), sp.diff(eq2(q1,q2,q3,q4),q3), "numpy")
    D8= sp.lambdify((q1,q2,q3,q4), sp.diff(eq2(q1,q2,q3,q4),q4), "numpy")
    
    return D1, D2, D3, D4, D5, D6, D7, D8


D1, D2, D3, D4, D5, D6, D7, D8 = differentiate(q1,q2,q3,q4)

q1_val = 105*pi/180
q2_val = 75*pi/180
q3_val = 0
q4_val = 0
DD1 = D1(q1_val,q2_val,q3_val,q4_val)
DD2 = D2(q1_val,q2_val,q3_val,q4_val)
DD3 = D3(q1_val,q2_val,q3_val,q4_val)
DD4 = D4(q1_val,q2_val,q3_val,q4_val)
DD5 = D5(q1_val,q2_val,q3_val,q4_val)
DD6 = D6(q1_val,q2_val,q3_val,q4_val)
DD7 = D7(q1_val,q2_val,q3_val,q4_val)
DD8 = D8(q1_val,q2_val,q3_val,q4_val)



A = np.array([[0,0,1,0],
        [0,     0,      0,  1],
        [DD1,  DD2,   DD3,  DD4],
        [DD5,  DD6,   DD7,  DD8]])

print(A)

print('trace:', np.trace(A))
print('eigenvals:', np.linalg.eigvals(A))

# print('\n22-23')
# def linearized_system(A,X0,t_range):
#     dt = t_range[1] - t_range[0]
#     sol = np.matmul(A,X0)
#     d1 = sol[2]
#     d2 = sol[3]
#     Xn = [X0[0]+X0[0]*d1*dt, X0[1]+X0[1]*d2*dt,d1,d2]
#     return Xn
# 
# 
# X0 = np.array([.35,.35,.35,.35])
# t_range = np.linspace(0,50,1000)
# 
# X = [X0]
# 
# for t in t_range:
#     #Xn = A*X[-1]
#     #Xn = np.matmul(A,X[-1])
#     Xn = linearized_system(A,X[-1],t_range)
#     X.append(Xn)
# 
# q1_f = X[-1][0]
# q2_f = X[-1][1]
# print(X)
# print(q1_f, q2_f)
# #sol = solve_ivp(linearized_system,[t_range[0],t_range[-1]],X0, method='Radau', t_eval = t_range)
# #print(np.shape(sol))

# print('problem 22,23')
# 
# D1, D2, D3, D4, D5, D6, D7, D8 = differentiate(q1,q2,q3,q4)
# 
# 
# X0 = np.array([.35,.35,.35,.35])
# X = [X0]
# t_range = np.linspace(0,50,10)
# for t in t_range:
#     q1_val = X[-1][0]
#     q2_val = X[-1][1]
#     q3_val = X[-1][2]
#     q4_val = X[-1][3]
# 
#     DD1 = D1(q1_val,q2_val,q3_val,q4_val)
#     DD2 = D2(q1_val,q2_val,q3_val,q4_val)
#     DD3 = D3(q1_val,q2_val,q3_val,q4_val)
#     DD4 = D4(q1_val,q2_val,q3_val,q4_val)
#     DD5 = D5(q1_val,q2_val,q3_val,q4_val)
#     DD6 = D6(q1_val,q2_val,q3_val,q4_val)
#     DD7 = D7(q1_val,q2_val,q3_val,q4_val)
#     DD8 = D8(q1_val,q2_val,q3_val,q4_val)
# 
#     A = np.array([[0,0,1,0],
#         [0,     0,      0,  1],
#         [DD1,  DD2,   DD3,  DD4],
#         [DD5,  DD6,   DD7,  DD8]])
#     Xn =  np.matmul(A,X[-1])
#     X.append(Xn)

print('problem 22,23')

def solve_acc(qs,dt):
    q1,q2,q3,q4 = qs
    shit1 = c5*q4*(2*q3+q4)*sin(q2)+ c6*sin(q1) + c7*sin(q1+q2) + 2*q3
    shit2 = d4*q3**2*sin(q2) + d5*sin(q1+q2) + 3*q4
    denom1 = c1+c2*cos(q2) - (c3+c4*cos(q2))*(d1+d2*cos(q2))/d3
    denom2 = c3+c4*cos(q2) - d3*(c1+c2*cos(q2))/(d1+d2*cos(q2))
    eq1 = (10*F*cos(q1)+ shit1 + shit2*(c3+c4*cos(q2))/d3)/denom1
    eq2 = (10*F*cos(q1) + shit1 + shit2*(c1+c2*cos(q2))/(d1+d2*cos(q2)))/denom2
    v1 = q3+eq1*dt
    v2 = q4+eq2*dt
    x1 = q1+v1*dt
    x2 = q2+v2*dt
    return x1,x2,v1,v2


X0 = np.array([.35,.35,.35,.35])
X = [X0]
t_range = np.linspace(0,50,1500000)
dt = t_range[1] - t_range[0]
for t in t_range:
    Xn = solve_acc(X[-1],dt)
    X.append(Xn)
    
print(X[-1])
import matplotlib.pyplot as plt
xs = [x[0] for x in X]
ys = [x[1] for x in X]
xs = xs[:-1]
ys = ys[:-1]
print(len(xs),len(ys),len(t_range))
plt.plot(t_range,xs)
plt.plot(t_range,ys)
plt.show()
