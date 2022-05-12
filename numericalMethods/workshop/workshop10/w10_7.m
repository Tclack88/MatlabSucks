% Landau equation (modeling small perturbations in a fluid)

dt = .01;
a = 0;
b = 10;
t_range = linspace(a,b,(b-a +dt)/dt);
A0 = .1;  %initial value
e = 1; % some initial epsilon guess
% Explicit Euler method
As1 = explicit_euler(@f,t_range, A0, e);

% 2nd order taylor - requires calculation of df/dt and df/dx (by hand)
As2 = taylor_2nd_order(@f, @dfdt, @dfdA, t_range,A0,e);

plot(t_range,As1,t_range,As2)
title('landau equation plotted plotted numerically')
xlabel('time')
ylabel('Amplitude')
legend('explicit euler','2nd order taylor')

% Let's try to vary e now
e_range = linspace(0,10,10);
final_A1 = []; %list of the final Euler A's (steady state amplitude A*)
final_A2 = []; %list of the final Taylor A's (steady state amplitude A*)
figure
hold on
for e = e_range;
	As1 = explicit_euler(@f,t_range, A0, e);
	As2 = taylor_2nd_order(@f, @dfdt, @dfdA, t_range,A0,e);
	final_A1(end+1) = As1(end);
	final_A2(end+1) = As2(end);
	plot(t_range,As1,t_range,As2)
end
hold off
xlabel('time')
ylabel('Amplitude')

% Now let's compare A* (steady state amplitude) against the choice of e
figure
hold on
plot(e_range,final_A1,'r+')


%%%%%%%%% functions %%%%%%%%%%%

function ret = f(A,e)
	t = 1; % tau, not time t
	g = 1;
	ret = e/t*A - g/t*A^3;
end

function ret = dfdA(A,e)
	% partial deriv of f wrt N
	t = 1; % tau, not time t
	g = 1;
	ret = e/t -3*g/t*A^2;
end

function ret = dfdt(A,e)
	% partial deriv of f wrt t
	ret = 0;
end

function As = explicit_euler(f,t_range,A0,e)
	As = zeros(size(t_range));
	As(1) = A0;
	dt = t_range(2)-t_range(1);
	for i = 1:length(t_range)-1
		As(i+1) = As(i) +  dt*f(As(i),e);
	end
end

function As = taylor_2nd_order(f, dfdt, dfdx,t_range,A0,e)
	As = zeros(size(t_range));
	As(1) = A0;
	dt = t_range(2)-t_range(1);
	for i = 1:length(t_range)-1;
		a = As(i);
		As(i+1) = a +  dt*f(a,e) + dt^2/2*(dfdt(a,e) + dfdx(a,e)*a);
	end
end
