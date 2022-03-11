% plot of force as a function of angle the force is applied
% for pulling a box along a ground with friction
u = .58;
g = 10;
m = 6;

F = @(t) u*m*g/(cos(t) +u*sin(t));
% We must look at the derivative of F, and find when the slope
% changes from being - to + 
dF = @(t) u*m*g*(sin(t) - u*cos(t))/(cos(t) + u*sin(t))^2

% set initial guesses
xl = .2;
xu = .6;

% I choose to use the span of the upper and lower bound as my criteria
% for convergence. Set to initially a large number.
gap = 10; 
n = 0;
while gap > .0000001
	xr = (xu + xl)/2;
	U = dF(xu);
	L = dF(xl);
	R = dF(xr);
	if R*L > 0
		xl = xr;
	elseif R*U > 0
		xu = xr;
	end
	n = n + 1;
	gap = xu - xl
	xr
end

fz = fzero(dF,.5);

sprintf('after %d cycles, a value of %f has been determined numerically as the lowest angle. \n the value determined by fzero is %f', n, xr, fz)

