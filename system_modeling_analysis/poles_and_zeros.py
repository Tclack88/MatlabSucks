import sympy as sym

x,s = sym.symbols('x,s')


def second_order_sys(K,Z,w,s):
    return K*s**2/(s**2+2*Z*w*s+w**2)

K,Z,w = 1,-1,4
system = second_order_sys(K,Z,w,s)

numer = sym.numer(system)
denom = sym.denom(system)

zeros = sym.solve(numer)
poles = sym.solve(denom)

print('system:')
print(system)
print('\nzeros:')
for z in zeros:
    print(z)
print('poles:')
for p in poles:
    print(p)
