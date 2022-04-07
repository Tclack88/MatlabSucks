% The Witch of Agnesi

x1 = [-1 0 1];
x2 = -1:.25:1;
xs = linspace(-1,1); % x doman for plotting

% Verify function works by checking the output of the above vals
f(x1)
f(x2)

hold on
plot(xs, f(xs))
lagrange_plot(x1, f(x1), xs)
lagrange_plot(x2, f(x2), xs)
legend('original', '3-point lagrange', '9-point lagrange')
title({'comparison of lagrange interpolation','colors don''t match because matlab sucks'})
hold off

% It appears that as more points are added, the lagrange interpolation
% polynomial gets "wavier" and does a poor job at approximating nearby points
% For many points, the interpolating polynomial is almost locally perpenddicular
% to the nodes leading to terrible local approximations 

x3 = [-1.0000, -0.8998, -0.6772, -0.3631, 0.0000, 0.3631, 0.6772, 0.8998, 1.0000]
figure
hold on
plot(xs, f(xs))
lagrange_plot(x3, f(x3), xs)
legend('original', 'lagrange interpolation', 'quadratic spline')
title('Approximation with Gauss–Lobatto–Legendre (GLL) nodes')
hold off

% These are apparently the more optimal choice for interpolation, it doesn't
% however seem like a great fiit to me. Local approximations are going to be
% off because of the near-perpendicular slope when the function is flat

% Compare to quadratic spline interpolation
figure
hold on
plot(xs, f(xs))
spline_plot(xs,x1,f(x1))
spline_plot(xs,x2,f(x2))
spline_plot(xs,x3,f(x3))
legend('original', 'quadratic spline-3 points', 'quadratic spline-9 points', 'quadratic spline-GL modes')
title({'comparison of spline interpolation','colors don''t match because matlab sucks'})
hold off

% Lagrange points with Gauss-Lobatto legendre modes seems to be the best
% Even compared to a range of spline techniques


%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%

function ret =  f(x)
       ret = 1./(1+25*x.^2);
end

function  lagrange_plot(nodes, y, x)
        L_f = zeros(length(x)); % L final, the result of summing the polynomials
        for i = 1:length(nodes);
                L = lagrange_L(i, nodes, x);
                L_f = L_f + L*y(i);
        end
        plot(x,L_f)
end



function [a,b,c] = quadratic_spline(xn,yn)
        % returns the coefficients for quadratic splines
	a = yn; % a values for quadratic spline equals the function height
	c = zeros(1,length(a)-1); % Initialize c coeffs
	c(1) = 0; % guess to get our final unknown; not necessary as 0
	h = xn(2:end) - xn(1:end-1); % h defined as the intervals
	for i = 2:length(xn)-1;
        	c(i) = (1/h(i))*(a(i+1)/h(i)-a(i)*((1/h(i-1))+(1/h(i)))+a(i-1)/h(i-1)-c(i-1)*h(i-1)); % derived in lecture, grabbed from workshop sample
	end
	b = zeros(1,length(xn)); % initialize b coeffs
	for i = 1:length(xn)-1;
        	b(i) = (a(i+1)-a(i))/h(i) -c(i)*h(i);
	end
end


function spline_plot(xs,xn, yn)
	[a,b,c] = quadratic_spline(xn,yn);
        f = zeros(length(xs)); % initialize the function to be plotted
        for x = 1:length(xs);
                for n = 2:length(xn);
                        if xn(n) >= xs(x); % if we're in the appropriate spline
                                f(x) = a(n-1) +b(n-1)*(xs(x)-xn(n-1))+c(n-1)*(xs(x)-xn(n-1))^2;
                                break
                        end
                end
        end
        plot(xs, f)
end

