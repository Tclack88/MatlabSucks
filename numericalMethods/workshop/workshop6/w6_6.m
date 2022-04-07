% Temperature of a lake and water depth
% z = depth (m), T = temperature (C)
z = [0 2.3 4.9 9.1 13.7 18.3 22.9 27.2];
T = [22.8 22.8 22.8 20.6 13.9 11.7 11.1 11.1];


% Quadratic Spline Interpolation
[a, b, c] = quadratic_spline(z,T);
xs = linspace(0,27.2);

% Lagrange interpolation
hold on
spline_plot(a,b,c,xs,z, T)
lagrange_plot(z, T, xs)
ms = spline(z,T,xs)
plot(xs, ms)
legend('original points', 'quadratic spline', 'lagrange', 'matlab built-in spline')
title({'comparison of interpolation techniques', 'colors don''t match up because matlab sucks'})
hold off


% The cubic spline has very little "waviness about the points
% This seems to me like a most realistic answer


function spline_plot(a,b,c,xs,zs, T)
	f = zeros(length(xs)); % initialize the function to be plotted
	for x = 1:length(xs);
		for z = 2:length(zs);
			if zs(z) >= xs(x); % if we're in the appropriate spline
				f(x) = a(z-1) +b(z-1)*(xs(x)-zs(z-1))+c(z-1)*(xs(x)-zs(z-1))^2;
				break
			end
		end
	end
	plot(zs, T, 'ro', xs, f)
	%legend('spline', 'original')
end

function [a,b,c] = quadratic_spline(z,T)
	% returns the coefficients for quadratic splines
a = T; % a values for quadratic spline equals the function height
c = zeros(1,length(a)-1); % Initialize c coeffs
c(1) = 0; % guess to get our final unknown; not necessary as 0
h = z(2:end) - z(1:end-1); % h defined as the intervals
for i = 2:length(z)-1;
	c(i) = (1/h(i))*(a(i+1)/h(i)-a(i)*((1/h(i-1))+(1/h(i)))+a(i-1)/h(i-1)-c(i-1)*h(i-1)); % derived in lecture, grabbed from workshop sample
end
b = zeros(1,length(z)); % initialize b coeffs
for i = 1:length(z)-1;
	b(i) = (a(i+1)-a(i))/h(i) -c(i)*h(i);
end
end

function  lagrange_plot(nodes, y, x)
	L_f = zeros(length(x)); % L final, the result of summing the polynomials
	for i = 1:length(nodes);
		L = lagrange_L(i, nodes, x);
		L_f = L_f + L*y(i);
	end
	plot(x,L_f)
end

