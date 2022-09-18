g = 9.8;
l = 1;
w = 10;
b = 1;
m = 3;

tspan = linspace(0,10);
X0 = [1.42; 0.01]; %

[t,y] = ode45(@calculate1, tspan,X0);
[t1,y1] = ode45(@calculate2, tspan,X0);
plot(t,y(:,1))
hold on
plot(t1,y1(:,1))
%Y
%length(Y(1,1:end-1))
%plot(tspan,Y(1,1:end-1))
%hold on;
%plot(tspan,Y(2,1:end-1))


function ret = calculate1(tspan,x)
	g = 9.8;
	l = 1;
	w = 10;
	b = 1;
	m = 3;

	A = [0 				1; 
	-97.8   -3*b/(m*l^2)];
	%w^2*(9/2*(g/(l*w))^2 -1)    -3*b/(m*l^2)];
	ret = A*x;
end

function ret = calculate2(tspan,x)
	g = 9.8;
	l = 1;
	w = 10;
	b = 1;
	m = 3;

	x1 = x(1);
	x2 = x(2);
	x1_dot = x2
	x2_dot = w^2*sin(x1)*cos(x1) - (3*g/2/l)*sin(x1) - (3*b)/(m*l^2)*x2
	ret = [x1_dot;x2_dot]
	%w^2*(9/2*(g/(l*w))^2 -1)    -3*b/(m*l^2)];
end

function ret = calculate(X0,tspan)
	g = 9.8;
	l = 1;
	w = 10;
	b = 1;
	m = 3;

	dt = tspan(2)-tspan(1);
	A = [0 				1; 
	w^2*(9/2*(g/(l*w))^2 -1)    -3*b/(m*l^2)];
	ret = zeros(2,length(tspan));
	ret(:,1) = X0 ;
	for i = 1:length(tspan);
		deriv = A*ret(:,1);
		ret(:,i+1) = ret(:,i) + deriv*dt;
	end
end
