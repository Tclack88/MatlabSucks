from math import exp, pi, sqrt, atan, log
w,Z = 2, .4
print('overshoot:')
print(exp(-Z*pi/sqrt(1-Z**2)))

print('peak time:')
print(pi/w)

print("rise time (10-90%):")
print((1/w)*(pi - atan(sqrt(1-Z**2)/Z)))

d = .02 # settling time, .02 is 2%
print(f'settling time ({round(d*100)}%):')
print(-log(d*sqrt(1-Z**2))/(Z*w))
