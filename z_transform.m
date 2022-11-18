%200*(s + 5)/((s + 0.1)*(s**2 + 2500)*(s**2 + 6*s + 10))
%g1 = zpk([-5],[-.1,50j,-50j,-3-j,-3+j],200)
% Tustin's method
g1 = zpk([],[0,-9],20)
c2d(g1,1/20,'tustin') % sampled at 20 hz (because T = 1/f)
