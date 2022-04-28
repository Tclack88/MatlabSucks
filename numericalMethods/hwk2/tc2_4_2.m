% Interpolation for a car accelerating - compare lagrange to cubic spline

t = [0 .5 1 2 3 8]; % time (s)
v = [0 11.1 13.2 14.8 16.2 16.1]; % velocity (m/s)


x_range = linspace(0,8); % domain for plotting


% cubic spline interpolation
figure
[a b c d] = cubic_spline(t,v);
spline_plot(a,b,c,d,x_range,t,v)
hold on

% lagrange polynomial interpolation
lagrange_plot(t, v, x_range)
title('Cubic spline plot vs Lagrange polynomial interpolation')
xlabel('time (s)')
ylabel('velocity (m/s)')
hold off
legend('known points','cubic spline', 'Lagrange polynomial', 'Location', 'southwest')


disp('Lagrange interpolating polynomials suffers from Runge''s phenomenom')
disp('(waviness between the points being interpolated.')
disp('With cubic spline, by ensuring the splines left and right of a point')
disp('share the same 2nd derivative, we fix this and get a better fit')


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



function [a,b,c,d] = cubic_spline(x, y)
	% returns the coefficients of the polynomial for the cubic spline
	% get a coeffs
	a = y(1:end-1); % From first assumption, a(i) = f(xi)

	h = x(2:end) - x(1:end-1); % h defined as intervals. Used for b,c,d
	% get c coeffs
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
	c = A\C; % c coefficients are the solution to this matrix equation

	% get b coeffs
	for i = 1:length(x)-2;
		b(i) =  1/h(i) *(a(i+1) - a(i)) - h(i)/3 *(2*c(i) + c(i+1));
	end

	% get d coeffs
	for i=1:length(x)-2;
		d(i) = (c(i+1) - c(i))/(3*h(i));
	end
	b(end+1) = 0; % fill inn with zeros to make the same length
	d(end+1) = 0; % fill in
end






function spline_plot(a,b,c,d,xs,x,y)
	% plot the cubic spline function from the coeffecients obtained
	f_c = zeros(1,length(xs)); %initialize funtion to be plotted
	for i = 1:length(xs)
		for j = 1:length(x)-1
			if x(j) > xs(i)
				j = j-1;
				break
			end
		end
		f_c(i) = a(j) + b(j)*(xs(i) - x(j)) + c(j)*(xs(i) - x(j))^2 + d(j)*(xs(i)-x(j))^3;
	end

	plot(x,y,'ro',xs,f_c)

end
