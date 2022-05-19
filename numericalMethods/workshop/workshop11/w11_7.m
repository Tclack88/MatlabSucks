%  IVP
% dx/dt = x^(1/5), with x(0)=0
% we can solve this by separation of variables:
% x^(-1/5) dx = dt
% integrating:
% 5/4 x^(4/5) = t + C
% applying BC, when t =0, x=0, so C= 0. Rearrange for x(t) yields:
% x(t) = (4t/5)^(5/4)

%Plot exact solution:
dt = .01; %step time
t = [0:dt:10];
x = f(t);
plot(t,x)

% Explicit euler
xs = explicit_euler(@df,t,dt);
hold on
plot(t,xs)
title('compare exact solution to explicit euler method')
legend('exact','explicit euler')
hold off

disp('explicit euler fails because x never grows. we feed it the initial condition x=0 and the next time step adds 0 to 0 and so on')

% Implicit euler
% To do implicit euler, we need to evaluate at f(x_n+1) instead of f(x_n)
% so xn+1 = xn + dt * (xn)^(1/5) becomes
%  xn+1 = xn + dt * (xn+1)^(1/5)
% so we need to use a foot finding with:
% h(p) = p - dt*p^(1/5) - xn    (where p = xn+1)

% Try Newton Raphson
% recall this is:  pn+1 = pn - h(p)/h'(p)
xs2 = zeros(size(t));
for n = 1:length(xs2)-1
	guess = 1; %probably bad, but it will always increase because it exp.
	xs2(n+1) = NR(guess,xs2(n),dt);
end
figure
plot(t,x,t,xs2);
title('compare exact to implicit euler')
legend('exact','implicit euler')

% Try fixed point iteration

xs3 = zeros(size(t));
for n = 1:length(xs2)-1
	guess = 1; %probably bad, but it will always increase because it exp.
	xs3(n+1) = fixed_point_iter(@fp1,guess,xs3(n),dt);
end
figure
plot(t,x,t,xs3);
title('compare exact to fixed point iteration')
legend('exact','fixed-point iteration')



%%%% functions %%%%%
function x = f(t)
	x = (4*t/5).^(5/4);
end

function ret = df(x)
	ret = x^(1/5);
end

function xs = explicit_euler(df,ts,dt)
	xs = zeros(size(ts));
	for n = 1:length(xs)-1
		xs(n+1) = xs(n) + dt* df(xs(n));
	end
end

function ret = NR(p,x0,dt)
	% Take guess value p and previous value x0 and try to get next val
	% iteratively with newton raphson method
	h  = 10;
	while h>.0001; % arbitrary choice
		h = p - dt*p^(1/5) -x0;
		dh = 1 - dt/5*p^(-4/5);
		p = p - h/dh;
	end
	ret = p;
end

function ret = fixed_point_iter(f,p,x0,dt)
	% p - initial guess
	% x0 - previous value
	% f - function being iterated over
	% dt - small time interval
	l = 10; % arbitrary large val
	max_iter = 100;
	n = 0;
	while l > .001 && n < max_iter
		g = f(x0,p,dt); %g = xo + dt*p^(1/5);
		l = abs(g);
		n = n+1;
	end
	ret = g;
end

function g = fp1(x0,p,dt)
	% one possible funciton for fixed point iteration
	g = x0 + dt*p^(1/5);
end

