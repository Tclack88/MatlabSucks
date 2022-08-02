

t_range = linspace(0,100,1000);
y0 = [1; 0; 0];

[t,y] = ode45(@f, t, y0);

plot(t,y(:,1))

function ret = f(t_range,x)
	A = [0 1 0; 0 0 1; -2/3 -5/3 -2/3];
	ret = A*x;
end
