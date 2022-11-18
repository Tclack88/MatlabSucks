from scipy.optimize import fsolve
import numpy as np
# An attempt to solve for Reynold's number with turbulent flow numerically
# in turbulent flow, zeta cannot be isolated so must be solved numerically
# currently (18Nov2022), it doesn't work, work on it later

def zeta_in_disguise(x):
    global k, d, Re
    return x + 2 * np.log10(2.51 * x / Re + k / (d * 3.71))

k = .046
Re = 625e3
Re =6.94e5
d = 250
#x = 1 / np.sqrt(zeta)
x = fsolve(zeta_in_disguise, 0)
print(x)
#let's test, if x is really the solution to the equation
print(-2 * np.log10(2.51 * x / Re + k / d / 3.71))
