clc;
clear all;
dt = .01;
t = [0:dt:1];

x = zeros(size(t));
x(1) = 1;
for n=1:length(t)-1
        k1=f(t(n), x(n));
	k2=f(t(n) + dt/2, x(n)+dt*k1/2);
        k3=f(t(n)+dt,x(n)+dt*k2)';
        x(n+1)=x(n) + (k1+4*k2 + k3)*dt/6;
end

[t1,x1] = ode45(@f,t,1);

plot(t,x,t1,x1)
title('3rd order runga kutta plot')
xlabel('t')
ylabel('x')
legend('my rk3 function','solution using ode45')
function ret = f(t, x)
	ret = (2/x)*cos(x^3) + 5*exp(-x/2);
end
