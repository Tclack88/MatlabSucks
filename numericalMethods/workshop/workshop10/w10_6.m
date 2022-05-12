% Population growth with logistic equation

dt = .01;
a = 0;
b = 10;
t_range = linspace(a,b,(b-a +dt)/dt);
N0 = 10;  %initial value
% Explicit Euler method
Ns1 = explicit_euler(@f,t_range, N0);

% 2nd order taylor - requires calculation of df/dt and df/dx (by hand)
Ns2 = taylor_2nd_order(@f, @dfdt, @dfdN, t_range,N0);

plot(t_range,Ns1,t_range,Ns2)
title('logistic curve plotted plotted numerically')
xlabel('time')
ylabel('population')
legend('explicit euler','2nd order taylor')



%%%%%%%%% functions %%%%%%%%%%%

function ret = f(N)
	r = 2;
	K = 1000;
	ret = r*N*(1-N/K);
end

function ret = dfdN(N)
	% partial deriv of f wrt N
	r = 2;
	K = 1000;
	ret = r*(1-2*N/K);
end

function ret = dfdt(N)
	% partial deriv of f wrt t
	ret = 0;
end

function Ns = explicit_euler(f,t_range,N0)
	Ns = zeros(size(t_range));
	Ns(1) = N0;
	dt = t_range(2)-t_range(1);
	for i = 1:length(t_range)-1
		Ns(i+1) = Ns(i) +  dt*f(Ns(i));
	end
end

function Ns = taylor_2nd_order(f, dfdt, dfdx,t_range,N0)
	Ns = zeros(size(t_range));
	Ns(1) = N0;
	dt = t_range(2)-t_range(1);
	for i = 1:length(t_range)-1;
		n = Ns(i);
		Ns(i+1) = n +  dt*f(n) + dt^2/2*(dfdt(n) + dfdx(n)*n);
	end
end
