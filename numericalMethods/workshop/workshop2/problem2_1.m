% Newton-Raphson method on the function:
% f =  40/x (1-exp(-8x/5)) = 9
% We do this by examining where the hyperfunction g is 0
% g = 40/x (1-exp(-8x/5)) - 9
% grpime = (40/x^2)*((8*x/5 + 1)*exp(-8*x/5) - 1)
g = 40/x (1-exp(-8*x/5)) - 9
gprime = (40/x^2)*((8*x/5 + 1)*exp(-8*x/5) - 1)



