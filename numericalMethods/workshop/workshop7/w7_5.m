% Integral approximation using Trapezoidal and Simpson's rule
% Over a single interval

% In order to calculate the error, we must compare it to the original
% integral. The actual function integrated from x=0 to 1 is given by:
% This can be done by hand/mentally but here is how to do it with Matlab

clc

actual = integral(@(x) f(x), 0,1)


Ts = [];
Ss = [];
for n=2:10
 	T_area = trapezoidal(@f, 0, 1, n);
 	S_area = simpsons(@f, 0, 1, n);
 	Ts(end+1) = T_area;
 	Ss(end+1) = S_area;
end


error_T = abs(actual-Ts)
error_S = abs(actual-Ss)

x_range = length(error_T);
% TODO fix this plott error BS from matlab
%plot(x_range, error_T, x_range, error_S)

%%%%%%%% functions %%%%%%%%%%%%

function ret = f(x)
	ret = 6*x.^3 + x.^2 -11*x + 4;
end

function area = trapezoidal(f, a, b, n)
	% Give a function f, lower bound a, upper bound b and 
	%number of trapezoids n, this function will approximate 
	%the area of a function using the trapezoidal rule
	range = linspace(a,b,n);
	delta = range(2) - range(1);
	mp = range(2:end-1); %mid points
	area = delta/2*(f(a) + 2*sum(f(mp)) + f(b));
end


function area = simpsons (f, a, b, n)
	% approximates area under curve using simplson's rule
	range = linspace(a,b,n);
	delta = range(2) - range(1);
	mp = range(2:end-1); %mid points
	mp_odd = mp(1:2:end);
	mp_even = mp(2:2:end);
	area = delta/3*(f(a) + 4*sum(f(mp_odd)) + 2*sum(f(mp_even)) + f(b));
end
