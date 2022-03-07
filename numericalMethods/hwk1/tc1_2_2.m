% use fixed point iteration to find roots where the function
% sin(2x) + (x^2)/3 -4 = -1
% so alternatively we are finding the function f(x)
% where f(x) = sin(2x) + (x^2)/3 -3 = 0
% with fixed point iteration, we approximate the next term from the previous
% x_i+1 = sqrt(3(3 - sin(2x_i)))
fx = @(x) sin(pi/180*2*x) + x^2/3 - 3; %pi/180 to get redians
gx = @(x) sqrt(3*(3 - sin(pi/180*2*x)));

guess1 = -2.5;
guess2 = 3.1;

h = 10; % for height. Initially relatively large number
x = guess1;
n = 0;
while h >  .0001
    x = gx(x);
    actual = fx(x);
    h = abs(actual); 	% abs because if our guess can be + or -. We really
    n = n + 1;		% just care about the magnitude of the "height"
end

sprintf('A close approximation of x=%f achieved after %d iterations',x, n)
