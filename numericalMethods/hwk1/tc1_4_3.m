% Use Newton-Raphson to plot the shape H(x,y) = 35
% Newton-Raphson as I know it is used for a single variable, so my strategy
% here is to set a range of x variables and determine what y must be in order
% to make H(x,y) = 35. This is straigtforward. A reminder H is given as
% H(x,y) = (y+6)^2 + 4x^2 -x^2y = 35
% so we are interested in some hyperfunction g(x)
% g(x,y) = (y+6)^2 + 4x^2 -x^2y - 35 = 0
x_range = -2:.1:2; % Save computation. reduce the gradations of y
y_range = []; % empty array to be appended to
for x=x_range
	yi = -abs(x); % It's somewhere between -3 and 0 and negative
	y = Newton(@H, @dH, x, yi, .0001);
	y_range = [y_range y];
end	

plot(x_range, y_range)
xlabel('x');
ylabel('y');
xlim([-3 3]);
ylim([-3 3]);
title("calculated x,y location of sea level at height 35")

% The plot we see is consistent with the tools created by matlab's
% built-in tools

function ret = H(x,y) 
	ret =  (y+6)^2 + 4*x^2 - y*x^2 - 35;
end

function ret = dH(x,y)
	%ret = 2*(x+6) + 8*x -2*x*y;
	ret = 2*(y+6) -x^2;
end

function ret = Newton(H, dH, x, yi, lim)
	% @H, @dH - the function and its derivative
	% xi - initial guess for x
	% y - the stationary y value
	% lim - small height for convergence criteria
	l = 10; %large initial value
	n = 0;
	while l > lim
		yi = yi - H(x,yi)/dH(x, yi);
		l = abs(H(x,yi));
		n = n + 1;
		if n > 50
			ret = inf; % set a number that won't be plotted
			break
		end
	end
	ret = yi;
end
