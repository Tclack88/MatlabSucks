import sympy as sym
import mpmath as mp
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import ZerosPolesGain
import cmath
import control

def highlight(word):
    space = ' '*round(len(word)+1)
    center = '#'+space+word+space+'#'
    end = '#'*len(center)
    print(f'\n{end}\n{center}\n{end}')


showall=0 # show all print statements
showall=1

s = sym.Symbol('s')
I = sym.Symbol('I')
Y = sym.Symbol('Y')
V = sym.Symbol('V')
O1 = sym.Symbol('O1')
O2 = sym.Symbol('O2')

r=.150
L=.700
m=1
J1 = 2.5e-3
J2 = 5e-4
k=3
c=.05
cp=.3
Ra=100
La=0
kt=2
ke=.5

# J1*s**2*O1 + c*s*O1 + 3*k*r**2*O1 - k*r**2*O2 - 2*k*r*Y - kt*V = 0
# J2*s**2*O2 + c*s*O2 + 3*k*r**2*O2 - k*r**2*O1 = 0
# m*s**2*Y + cp*s*Y + 4*k*Y - 2*k*r*O1 - 2*k*r*O2 = 0
# La*s*I + Ra*I + ke*s*O1 = V

if showall:
    ##### PART 1 #####
    highlight('PART 1')
    solution = sym.solve((
        J1*s**2*O1 + c*s*O1 + 3*k*r**2*O1 - k*r**2*O2 - 2*k*r*Y - kt*I,
        J2*s**2*O2 + c*s*O2 + 3*k*r**2*O2 - k*r**2*O1 -2*k*r*Y,
        m*s**2*Y + cp*s*Y + 4*k*Y - 2*k*r*O1 - 2*k*r*O2,
        La*s*I + Ra*I + ke*s*O1 - V), (V,Y,I,O1,O2))

    print(sym.simplify(solution[Y]/solution[V]))
    # numerator s term:    720
    # denominator s term: 151308o/10 = 151308
    highlight('Q1')
    print('numerator s term:    720')
    highlight('Q2')
    print('# denominator s term: 151308o/10 = 151308')
    
    ## Q3
    solution = sym.simplify((.01*s**5 + 1.203*s**4 + 25.34*s**3 + 183.858*s**2 + 619.08*s+1383.48)/(s**5+124.3*s**4+2935.2*s**3+20173.8*s**2+67194*s+151308))
    print(solution,'\n')
    highlight('Q3')
    print('5th order max in numerator => 5')

    highlight('Q4, Q5')
    
    #### Errors with sympy invlap #####
    # t = sym.Symbol('t') # define a t symbol for laplace transform
    # #step_resp = sym.simplify(solution*(1/s))
    # def invL(F):
    #     return sym.inverse_laplace_transform(F,s,t)
    # 
    # #inv_lap = sym.inverse_laplace_transform(step_resp, s, t)
    # print(sym.simplify(solution*(1/s),rational=True))
    # print()
    # print(solution*(1/s))
    # step_resp = sym.simplify(solution*(1/s),rational=True)
    # inv_lap = invL(step_resp).evalf()
    # print(inv_lap)
    
    # workaround as sympy inv laplace is erroring
    step_resp = sym.simplify(solution*(1/s),rational=True)
    #print(step_resp)
    
    def f(s): # the step response (printed above for convenience) as a function for mpmath
        return (10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(100*s*(10*s**5 + 1243*s**4 + 29352*s**3 + 201738*s**2 + 671940*s + 1513080))
    
    t = np.linspace(0.001,12,100)
    
    G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
    glist = list(G)
    c = np.array(glist, dtype=float)
    print(c[0],c[-1]) #at t=0 and t=12
    #0.009960476716492804 0.009143468959732321
    plt.plot(t,c)
    plt.show()

    zeros = sym.solve(-(10*s**4 + 1003*s**3 + 4470*s**2 + 13215*s + 32400))
    poles = sym.solve(25*(10*s**5 + 1243*s**4 + 29352*s**3 + 201738*s**2 + 671940*s + 1513080))
    print('zeros')
    for z in zeros:
      print(z.evalf())
    print('poles')
    for p in poles:
      print(p.evalf())
    
    # two complex poles:
    # -1.44844812267064 - 3.40154098071737*I
    # -1.44844812267064 + 3.40154098071737*I
    # so real part is -1.45
    highlight('Q6')
    print(-1.448448)
    highlight('Q7')
    print(abs(-1.44844812267064 + 3.40154098071737j))
    # 3.697
    highlight('Q8')
    print('negative real and imag means stable and oscillatory')
    


    ##### PART 2 #####
    highlight('PART 2')
    
    #Gol was "solution" from part1 Q3. Gol = open loop transfer function
    Gol = sym.simplify((.01*s**5 + 1.203*s**4 + 25.34*s**3 + 183.858*s**2 + 619.08*s+1383.48)/(s**5+124.3*s**4+2935.2*s**3+20173.8*s**2+67194*s+151308))
    
    # Gcl - closed loop transfer function
    
    Ks =  [750,1000,2000,3000,4000]
    for K in Ks:
        fac = K/(s-5)
        Gcl = fac*Gol/(1+fac*Gol)
        print(K)
        print(sym.simplify(Gcl,rational=True))
    
    # 750
    # 15*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(20*s**6 + 2536*s**5 + 64319*s**4 + 490056*s**3 + 2084370*s**2 + 5592960*s + 5621400)
    # 1000
    # 10*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(10*s**6 + 1293*s**5 + 35167*s**4 + 308378*s**3 + 1501830*s**2 + 4344180*s + 6269400)
    # 2000
    # 20*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(10*s**6 + 1393*s**5 + 47197*s**4 + 561778*s**3 + 3340410*s**2 + 10534980*s + 20104200)
    # 3000
    # 30*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(10*s**6 + 1493*s**5 + 59227*s**4 + 815178*s**3 + 5178990*s**2 + 16725780*s + 33939000)
    # 4000
    # 40*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(10*s**6 + 1593*s**5 + 71257*s**4 + 1068578*s**3 + 7017570*s**2 + 22916580*s + 47773800)
    # 
    nums=[15*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480),
    10*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480),
    20*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480),
    30*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480),
    40*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)]
    
    denoms =[ 
    (20*s**6 + 2536*s**5 + 64319*s**4 + 490056*s**3 + 2084370*s**2 + 5592960*s + 5621400),
    (10*s**6 + 1293*s**5 + 35167*s**4 + 308378*s**3 + 1501830*s**2 + 4344180*s + 6269400),
    (10*s**6 + 1393*s**5 + 47197*s**4 + 561778*s**3 + 3340410*s**2 + 10534980*s + 20104200),
    (10*s**6 + 1493*s**5 + 59227*s**4 + 815178*s**3 + 5178990*s**2 + 16725780*s + 33939000),
    (10*s**6 + 1593*s**5 + 71257*s**4 + 1068578*s**3 + 7017570*s**2 + 22916580*s + 47773800)]
    
    highlight('Q9')
    print('6th order in denomical max => 6')
    
    highlight('Q10-12 : poles')
    for i,K in enumerate(Ks):
        print()
        print(K)
        poles = sym.solve(denoms[i])
        print('poles')
        for j,p in enumerate(poles):
            if j == 2 and K in [750,2000,4000]:
                print(f'\t{K}\t--------> {p.evalf()}')
            else:
                print(p.evalf())
    
    highlight('Q13-15: poles')
    for i,K in enumerate(Ks):
        print()
        print(K)
        zeros = sym.solve(nums[i])
        print('zeros')
        for j,z in enumerate(zeros):
            if j == 2 and K in [1000,2000,3000]:
                print(f'\t{K}\t--------> {z.evalf()}')
            else:
                print(z.evalf())
    
    highlight('Q16')
    
    K_4000 = 40*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(10*s**6 + 1593*s**5 + 71257*s**4 + 1068578*s**3 + 7017570*s**2 + 22916580*s + 47773800)
    # The below didn't work but it was a good try
    # coeff_vars = sym.var('A:E')
    # poles = sym.solve(sym.denom(K_4000))
    # poles_eval = [s - p.evalf() for p in poles]
    # # build equation of manually separating the partial fraction decomposition
    # eq = None
    # for i,n in enumerate(coeff_vars):
    #   if not eq:
    #     eq = n/poles_eval[i]
    #   else:
    #     eq += n/poles_eval[i]
    # abcde = sym.solve((sym.numer(K_4000), eq),(A,B,C,D,E))
    
    k4000_apart = sym.apart(K_4000,full=True).evalf()
    #print(k4000_apart)
    easy_denoms =  + 0.000371254389165188/(s + 95.7790570761686) + 35.1699224573921/(s + 42.7395905270505)
    print(easy_denoms)
    denom_with_8p73 = (-11.2464560034781 + 503048965888709690038037091641218893842737672*(-8.73239211709172 - 2.23520508415252j)**3/9442704057835804252257989726393814888969168617 - 4.13726122181375j + 87328480679641149818341306193204948154960*(-8.73239211709172 - 2.23520508415252j)**5/9442704057835804252257989726393814888969168617 + 770276017899818288604411067806486052450184*(-8.73239211709172 - 2.23520508415252j)**4/555453179872694367779881748611400875821715801 + 4963564646332212836733215351596015097695522128*(-8.73239211709172 - 2.23520508415252j)**2/9442704057835804252257989726393814888969168617)/(s + 8.73239211709172 + 2.23520508415252j) + (-11.2464560034781 + 4963564646332212836733215351596015097695522128*(-8.73239211709172 + 2.23520508415252j)**2/9442704057835804252257989726393814888969168617 + 770276017899818288604411067806486052450184*(-8.73239211709172 + 2.23520508415252j)**4/555453179872694367779881748611400875821715801 + 87328480679641149818341306193204948154960*(-8.73239211709172 + 2.23520508415252j)**5/9442704057835804252257989726393814888969168617 + 4.13726122181375j + 503048965888709690038037091641218893842737672*(-8.73239211709172 + 2.23520508415252j)**3/9442704057835804252257989726393814888969168617)/(s + 8.73239211709172 - 2.23520508415252j)
    print(denom_with_8p73)
    
    denom_with_1p65 = (1.84739075935892 - 6.30782192930868j + 770276017899818288604411067806486052450184*(-1.65828408129873 - 3.40787658559749j)**4/555453179872694367779881748611400875821715801 + 87328480679641149818341306193204948154960*(-1.65828408129873 - 3.40787658559749j)**5/9442704057835804252257989726393814888969168617 + 503048965888709690038037091641218893842737672*(-1.65828408129873 - 3.40787658559749j)**3/9442704057835804252257989726393814888969168617 + 4963564646332212836733215351596015097695522128*(-1.65828408129873 - 3.40787658559749j)**2/9442704057835804252257989726393814888969168617)/(s + 1.65828408129873 + 3.40787658559749j) + (1.84739075935892 + 4963564646332212836733215351596015097695522128*(-1.65828408129873 + 3.40787658559749j)**2/9442704057835804252257989726393814888969168617 + 503048965888709690038037091641218893842737672*(-1.65828408129873 + 3.40787658559749j)**3/9442704057835804252257989726393814888969168617 + 87328480679641149818341306193204948154960*(-1.65828408129873 + 3.40787658559749j)**5/9442704057835804252257989726393814888969168617 + 770276017899818288604411067806486052450184*(-1.65828408129873 + 3.40787658559749j)**4/555453179872694367779881748611400875821715801 + 6.30782192930868j)/(s + 1.65828408129873 - 3.40787658559749j)
    print(denom_with_1p65)
    
    highlight('Q17')
    print('stable and oscillatory, because negative real and imaginaty')
    
    highlight('Q18-22')
    
    
    Ks =  [750,1000,2000,3000,4000]
    Gol = sym.simplify((.01*s**5 + 1.203*s**4 + 25.34*s**3 + 183.858*s**2 + 619.08*s+1383.48)/(s**5+124.3*s**4+2935.2*s**3+20173.8*s**2+67194*s+151308))
    for K in Ks:
        print(K)
        fac = K/(s-5)
        Gcl = fac*Gol/(1+fac*Gol)
        step_resp = sym.simplify(Gcl*(1/s))
    
        def f(s): # the step response (found from the print statement above for convenience) as a function for mpmath
            step_dict={750: 15*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(s*(20*s**6 + 2536*s**5 + 64319*s**4 + 490056*s**3 + 2084370*s**2 + 5592960*s + 5621400)),
            1000:10*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(s*(10*s**6 + 1293*s**5 + 35167*s**4 + 308378*s**3 + 1501830*s**2 + 4344180*s + 6269400)),
            2000:20*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(s*(10*s**6 + 1393*s**5 + 47197*s**4 + 561778*s**3 + 3340410*s**2 + 10534980*s + 20104200)),
            3000:30*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(s*(10*s**6 + 1493*s**5 + 59227*s**4 + 815178*s**3 + 5178990*s**2 + 16725780*s + 33939000)),
            4000:40*(10*s**5 + 1203*s**4 + 25340*s**3 + 183858*s**2 + 619080*s + 1383480)/(s*(10*s**6 + 1593*s**5 + 71257*s**4 + 1068578*s**3 + 7017570*s**2 + 22916580*s + 47773800))}
            return step_dict[K]
        
        t = np.linspace(.1,10,100)
        G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
        glist = list(G)
        c = np.array(glist, dtype=float)
        print('\t i-r=',c[-1] -1)

highlight('Part 3')

highlight('Q23')

K = 4.1
w = 78.4
Z = 1.9

def f(s):
    # return the step response as a function for mpmath
    return (K*w**2/(s**2 + 2*Z*w*s + w**2))*1/s

t = np.linspace(.1,10,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
glist = list(G)
c = np.array(glist, dtype=float)
print(max(c))
plt.plot(t,c)
plt.title('Q23 - max value ~ 4.1')
plt.show()


#import sys; sys.exit(0)

highlight('Q24')
Ve = 4*s**4 + 401.2*s**3 + 1788*s**2 + 5286*s + 12960
V = s**5 + 124.3*s**4 + 2935.2*s**3 + 20173.8*s**2 + 67194*s + 151308
Tf = sym.simplify(Ve/V)
print(Tf)

zeros = sym.solve(Ve)
print('zeros')
for z in zeros:
    print(z)
print('poles')
poles = sym.solve(V)
for p in poles:
    print(p)

print(sym.simplify((s-poles[-1])*(s-poles[-2])))
# last value is 13.6684830075678 = w^2
wn = 13.6684830075678**.5
print(f'wn = {wn}')
print(f'K = {1/13.6684830075678}')
print(f'zeta = {2.89689624534129/(wn*2)}') 

def f(s):
    wn = 3.6970911548902605
    K = 0.07316100839034823
    zeta = 0.39178047334719796

    # return the step response as a function for mpmath
    return (K*w**2/(s**2 + 2*Z*w*s + w**2))*1/s

t = np.linspace(.1,10,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
glist = list(G)
c = np.array(glist, dtype=float)
print(max(c))
plt.plot(t,c)
plt.show()


import sys; sys.exit(0)
# plot to check normal (unsimplified) step response
# def f(s):
#     Ve = 4*s**4 + 401.2*s**3 + 1788*s**2 + 5286*s + 12960
#     V = s**5 + 124.3*s**4 + 2935.2*s**3 + 20173.8*s**2 + 67194*s + 151308
#     return (Ve/V)/s
# 
# t = np.linspace(.1,10,100)
# G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
# glist = list(G)
# c = np.array(glist, dtype=float)
# print(max(c))
# plt.plot(t,c)
# plt.show()


def f(s):
    #r1 = s - (-1.44844812267064 - 3.40154098071737j)
    #r2 = s - (-1.44844812267064 + 3.40154098071737j)
    r1 = s - poles[-1]
    r2 = s = poles[-2]
    #r3 = s - (-5.84265829258361)
    #return 1/r3
    return (12960/(r1*r2))/s

t = np.linspace(.1,10,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t)
glist = list(G)
c = np.array(glist, dtype=float)
print(max(c))
plt.plot(t,c)
plt.show()

# num = [4,401.2,1788,5286,12960]
# denom = [1,124.3,2935.2,20173.8,67194,151308]
# F = control.tf(num, denom)
# control.rlocus(F)
# plt.grid()
# plt.show()

