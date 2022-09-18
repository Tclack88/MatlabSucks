clear
close all

m = 1;
g = 9.8;
L = 1;
b = 1;

%% Derivative function
odefn = @(t,x) [x(2);
                (-b*x(2)-m*g*L*sin(x(1)))/(m*L^2)];

%% Initial state and timespan
X0 = [deg2rad(90),0];
timespan = [0,10];

%% Simulation
[t,X] = ode45(odefn,timespan,X0);

%% Plot x from theta
figure()
theta = X(:,1);
x = L*sin(theta);
plot(t,x)