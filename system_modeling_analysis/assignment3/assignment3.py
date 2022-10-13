from scipy import signal
import matplotlib.pyplot as plt
import sympy as sym
import cmath
import numpy as np
import mpmath as mp

s = sym.Symbol('s')
t = sym.Symbol('t')


def highlight(word):
    space = ' '*round(len(word)+1)
    center = '#'+space+word+space+'#'
    end = '#'*len(center)
    print(f'\n{end}\n{center}\n{end}')

showall = 0

if showall:
    highlight('Part 1')
    highlight('Q1-2')
    
    num = [7.2, 720, 3888]
    denom = [1, 124.3, 2935.2, 20173.8, 67194, 151308, 0]
    
    sys = signal.TransferFunction(num,denom)
    w,mag,phase = signal.bode(sys)
    plt.figure()
    plt.semilogx(w,mag)
    plt.grid()
    plt.title('mag plot')
    
    plt.figure()
    plt.semilogx(w,phase)
    plt.title('phase plot')
    plt.grid()
    plt.show()
    
    print('low gain at high frequencies')
    print('low (negative) phase at higher frequencies')
    
highlight('Q3')

Gol = (7.2*s**2 + 720*s + 3888)/(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s)

lap_sin = sym.laplace_transform(sym.sin(t), t, s, noconds=True)
#print(lap_sin)
A = 500
def lap_f(s):
    return sym.simplify(A*lap_sin*Gol,rational=True)

#f = lap_f(s)
#print(f)

# Failed to do inverse laplace transform symbolically
# F = sym.inverse_laplace_transform(f,s,t)
# print(F)
# #sym.plot(f(s).
# import sys; sys.exit(0)
def f(s):
    # printed from above
    #(because there's an error if I don't write it explicitly)
    return 36000*(s**2 + 100*s + 540)/(s*(s**2 + 1)*(10*s**5 + 1243*s**4 + 29352*s**3 + 201738*s**2 + 671940*s + 1513080))
t_range = np.linspace(0.001,50,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
glist = list(G)
c = np.array(glist, dtype=float)
plt.plot(t_range,c)
plt.title('Q3')
plt.show()
print('steady state value', c[-1],'\nalso looks like it\'s off by (Q3) -90')


highlight('Q4')
A = 10
def lap_f(s):
    return A*lap_sin*Gol
#print(lap_f(s))

def f(s):
    # printed from above
    #(because there's an error if I don't write it explicitly)
    return 10*(7.2*s**2 + 720*s + 3888)/((s**2 + 1)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))


t_range = np.linspace(0.001,50,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
glist = list(G)
c = np.array(glist, dtype=float)
plt.plot(t_range,c)
plt.title('Q4')
plt.show()
print('steady state value (Q4)', c[-1],'\nalso looks like it\'s off by -90')


highlight('Q5')
A = 1000
w = .5
lap_sin = sym.laplace_transform(A*sym.sin(w*t), t, s, noconds=True)
def lap_f(s):
    return lap_sin*Gol
#print(lap_f(s))

def f(s):
    # printed from above
    #(because there's an error if I don't write it explicitly)
    return 2000.0*(7.2*s**2 + 720*s + 3888)/((4.0*s**2 + 1)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))

t_range = np.linspace(0.001,100,500)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
glist = list(G)
c = np.array(glist, dtype=float)
plt.plot(t_range,c)
plt.title('Q5')
plt.show()
print('steady state value', c[-1],'\nalso looks like it\'s off by -90')


highlight('Q6')
w = 50
lap_sin = sym.laplace_transform(A*sym.sin(w*t), t, s, noconds=True)
def lap_f(s):
    return lap_sin*Gol
#print(lap_f(s))

def f(s):
    # printed from above
    #(because there's an error if I don't write it explicitly)
    return 50000*(7.2*s**2 + 720*s + 3888)/((s**2 + 2500)*(s**6 + 124.3*s**5 + 2935.2*s**4 + 20173.8*s**3 + 67194*s**2 + 151308*s))

t_range = np.linspace(0.001,50,100)
G = map( lambda x: mp.invertlaplace(f, x, method = 'dehoog', dps = 10, degree = 18), t_range)
glist = list(G)
c = np.array(glist, dtype=float)
plt.plot(t_range,c)
plt.show()
plt.title('Q6')
max_val = max(c)
print(f'max deviation (max val - ss): {max_val} - {c[-1]} = {abs(max_val - c[-1])}')

highlight('Q7')

highlight('Part 2')
Gol = (.35*(s+5.7))/(s*(s+5.8)*(s+1.4+3.4j)*(1+1.4-3.4j))
w = 50
K = 15
