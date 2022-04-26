% Interpolation for a car accelerating

t = [0 .5 1 2 3 8]; % time (s)
v = [0 11.1 13.2 14.8 16.2 16.1]; % velocity (m/s)


x_range = linspace(0,8); % domain for plotting

%plot(t,v)

% cubic spline interpolation
% spline_y_range = spline(t,v,x_range)
% 
% plot(t,v,'ko','MarkerSize',5)
% hold on
% plot(x_range, spline_y_range)

% lagrange polynomial interpolation

lagrange_plot(t, v, x_range)
title('Cubic spline plot vs Lagrange polynomial interpolation')
xlabel('time (s)')
ylabel('velocity (m/s)')
hold off
%legend('known points','cubic spline', 'Lagrange polynomial', 'Location', 'southwest')
% 
figure
[a b c d] = cubic_spline(t,v);
hold on
spline_plot(a,b,c,d,x_range, t,v);
%[a1 b1 c1] = quadratic_spline(t,v)
%q_spline_plot(a1,b1,c1, x_range, t,v)
pp = csape(t,v,'variational');
%pp = spline(t,v); % not a natural cubic spline
pp.coefs

a1 = pp.coefs(1:end,4)
b1 = pp.coefs(1:end,3)
c1 = pp.coefs(1:end,2)
d1 = pp.coefs(1:end,1)

spline_plot(a1,b1,c1,d1,x_range,t,v)
hold off


%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%

function L = lagrange_L(i,nodes,x)
	% A provided script
        % i - ith lagrange polynomial to be calculated
        % nodes: known x vals of the original function
        % x_range domain (many points for interpolation)

L = ones(length(x),1);

for k = 1:length(x);
    for j = 1:length(nodes);
        if j~=i
            L(k) = L(k)*(x(k)-nodes(j))/(nodes(i)-nodes(j));
        end
    end
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


% function [a,b,c] = quadratic_spline(z,T)
%         % returns the coefficients for quadratic splines
% a = T; % a values for quadratic spline equals the function height
% c = zeros(1,length(a)-1); % Initialize c coeffs
% c(1) = 0; % guess to get our final unknown; not necessary as 0
% h = z(2:end) - z(1:end-1); % h defined as the intervals
% for i = 2:length(z)-1;
%         c(i) = (1/h(i))*(a(i+1)/h(i)-a(i)*((1/h(i-1))+(1/h(i)))+a(i-1)/h(i-1)-c(i-1)*h(i-1)); % derived in lecture, grabbed from workshop sample
% end
% b = zeros(1,length(z)); % initialize b coeffs
% for i = 1:length(z)-1;
%         b(i) = (a(i+1)-a(i))/h(i) -c(i)*h(i);
% end
% end


% function q_spline_plot(a,b,c,xs,zs, T)
%         f = zeros(length(xs)); % initialize the function to be plotted
%         for x = 1:length(xs);
%                 for z = 2:length(zs);
%                         if zs(z) >= xs(x); % if we're in the appropriate spline
%                                 f(x) = a(z-1) +b(z-1)*(xs(x)-zs(z-1))+c(z-1)*(xs(x)-zs(z-1))^2;
%                                 break
%                         end
%                 end
%         end
%         plot(zs, T, 'ro', xs, f)
%         %legend('spline', 'original')
% end
% 
% 

function [a,b,c,d] = cubic_spline(x, y)
	% returns the coefficients of the polynomial for the cubic spline
	a = y(1:end-1); % From first assumption, a(i) = f(xi)
	h = x(2:end) - x(1:end-1); % h defined as intervals. Used for b,c,d
	% c coefficients are the solutions of a certain tridiagonal matrix
	A = zeros(length(y)-1);
	A(1) = 1;
	A(end) = 1;
	for i = 2:length(A)-1;
		A(i,i-1) = h(i-1);
		A(i,i) = 2*(h(i-1)+h(i));
		A(i,i+1) = h(i);
	end
	C(1) = 0;
	C(end) = 0;
	C = zeros(length(A),1);
	for i = 2:length(C)-1;
		C(i) = (3/h(i))*(a(i+1)-a(i)) + (3/h(i-1))*(a(i-1) - a(i));
	end
	c = A\C % c coefficients are the solution to this matrix equation
	b = zeros(1,length(a)-1);
	for i = 1:length(b);
		b(i) = (1/h(i))*(a(i+1) - a(i)) - (h(i)/3)*(2*c(i) + c(i+1));
		%b = (1./h(1:end-1).*(a(2:end)-a(1:end-1)) - h(1:end-1)/3 .*(2*c(1:end-1) + c(2:end))
	end
	d = zeros(1,length(b));
	for i=1:length(d);
		d(i) = (c(i+1) - c(i))/(3*h(i));
	end
	a
	b
	c
	d
	%b(end) = 0
	%d(end) = 0
end

function spline_plot(a,b,d,c,xs, x, y)
	% plot the cubic spline function from the coeffecients obtained
	%f_c = zeros(size(xs)); %initialize funtion to be plotted
	for i = 1:length(xs);
		for j = 2:length(x)-1;
			if x(j) >= xs(i); % if we're in the appropriate spline
				f_c(i) = a(j-1) + b(j-1)*(xs(i) - y(j-1)) + c(j-1)*(xs(i) - y(j-1))^2 + d(j-1)*(xs(i) - y(j-1))^3;
				break
			end
		end
	end
	disp('debug')
	length(xs)
	length(f_c)
	%f
	plot(x,y,'ro',xs,f_c)
end
