tspan = linspace(0,90);
y0 = [0; 0; 0; 0; 0; 0];

[t,y] = ode45(@f, tspan, y0);

t(1)
t(end)
yplot = y(:,5)
y(1)
t(end)
plot(t, yplot)
title('y as a function of time')

function ret = f(t,x)
	V = 15*sin(.2*pi*t);
	ia = (2*V - x(2))/32;

	A = [0      1     0      0        0      0;
	    -1350   -2000 450    0        9000   0;
	     0      0     0      1        0      0;
	     300    0     -900   -1333.3  6000   0;
	     0      0     0      0        1      0;
	     1.385    0     1.385    0        -27.7  -.23];

    	c = [0;
	     10000*ia;
	     0;
	     0;
	     0;
	     0];
	ret = A*x + c;
end
