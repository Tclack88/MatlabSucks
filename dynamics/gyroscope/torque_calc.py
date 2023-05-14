# large disk side:
m1 = 1.747



# smal disk side:
m2 = .9 
m3 = .03
g = 9.8


small_CW = (48/13.2*8.1)/100
large_CW = (48/13.2*7)/100
large_wheel = (48/13.2*3.6)/100

print(large_wheel)
print(large_CW)
print(small_CW)

T2 = (m1*g*large_wheel)

T1 = (g*(m2*large_CW + m3*small_CW))

r1o = .254/2
r1i = .0093/2
L =  (T1-T2)/.1
w = L/(m1*(r1o-r1i)**2/2)
print(w)


