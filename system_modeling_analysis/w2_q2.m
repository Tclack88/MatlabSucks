tspan = linspace(0,.2);
y0 = [0; 0; 0; 0];

[t,y] = ode45(@f, tspan, y0);

plot(t, y(:,1))
hold on
plot(t, y(:,3))

function ret = f(t,x)
	b = .004;
	I1 = .001;
	I2 = .002;
	k = 10;
	amp = 10;
	freq = 100;
	M = amp*sin(2*pi*freq*t);
	A = [0       1     0    0;
	    -k/I1 -b/I1 k/I1 b/I1;
	     0       0     0    1;
	    k/I2  b/I2 -k/I2 -b/I2];
    	c = [0; M/I1; 0;0];
	ret = A*x + c;
end
