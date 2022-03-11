% out of curiosity, I will plot each successive terms
% against the number of new iterations to observe the behavior
fx = @(x) sqrt(3*(3-sin(2*x)))

xi = 3.1 % initial guess
xs = [];
n = [1:50];
for i = n
	xi = fx(xi);
	xs(end +1) = xi;
end

plot(n,xs)
disp(xs)
% the term oscillates around the value it's converging to
