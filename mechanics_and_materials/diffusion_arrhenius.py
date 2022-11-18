from math import exp, sqrt, log
from scipy.special import erfinv, erf

####### Diffusion #######
# Q: units ev/atam
# D: units m^2/s
# R: 8.314J/molK or 6.26e-5 ev/atomK
R1= 8.314 #J/molK
R3 =  6.26e-5 #ev/atomK # this is th WRONG value
R2 =  8.62e-5 #ev/atomK

# D = Do*exp(-Q/(R*T))
# (cx - co)/(cs - co) = 1 - erf(x/(2*sqrt(D*t)))   

def D_calc(Q,T,R,Do):
    D = Do*exp(-Q/(R*T))
    return D

def cx_calc(co,cs,x,D,t):
    #cx concentration at depth x after time t
    #co initial concentration throughout the material
    #cs concentration kept at the surface during the treatment
    cx = co + (cs-co)*(1- erf(x/(2*sqrt(D*t))))
    return cx

def t_calc(cs,cx,co,x,D):
    t = ((x/2) * 1/(erfinv(1 - (cx-co)/(cs-co))))**2/D
    return t

def x_calc(cs,cx,co,t,D):
    # x distance returned in m (assuming D is given in m^2/s)
    x = 2*sqrt(D*t)*erfinv(1 - (cx-co)/(cs-co))
    return x



t = 3600 # seconds
cs = 0
co = .12
cx = .06

T = 900+273 # Kelvin
#T = 400+273 # Kelvin
#T = 800+273 # Kelvin

x = .1e-3

#Carbon
Q = .83 # ev/atom
Do = 1e-6 # m^2/s

#Chromium
#Q=2.49
#Do=2e-4


D = D_calc(Q,T,R2,Do)
x = x_calc(cs,cx,co,t,D)
#print(x)


co = .2
x = 10e-6
T = 750+273
Do = 9.5e-6 #  9.5 mm^2/s -> m^2/s
Q = .159 #J/mol
D = D_calc(Q,T,R1,Do)
t = x**2/D
print(t)
#t_calc(cs,cx,co,x,D)


import sys; sys.exit(0)
# texbook example check

# x=.5e-3
# T = 950+273
# cx=.8
# cs = 1.2
# co = .25
# 
# D = 1.6e-11
# t = t_calc(cs,cx,co,x,D)
# print(t, t/3600)
# 
# 
# Q= 130000
# T = 550 +273
# Do = 1.2e-4
# D = D_calc(Q,T,R1,Do)
# print(D)



## Arrhenius stuff
def Q_calc(T1,T2,k1,k2,R):
    num = -R*(log(k1) - log(k2))
    denom = 1/T1 - 1/T2
    return num/denom

def k_calc(T1,T2,k2,Q,R):
    k1 = log(k2) - Q/R*(1/T1 - 1/T2)
    k1_2 = k2*exp(-Q/R*(1/T1 - 1/T2))
    return exp(k1)

T1 = 180+273
T2 = 140+273
T3 = 227+273
k1 = 30
k2 = 90


Q=Q_calc(T1,T2,k1,k2,R1)
print(Q)

k3 = k_calc(T3,T2,k2,Q,R1)
print(k3)


