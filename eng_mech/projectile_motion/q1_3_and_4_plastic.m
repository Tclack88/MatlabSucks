clc;
clear;

%%% Plastic %%%
% get data from fig
open('plastic-XY.fig');
a = get(gca,'Children');
x_real = get(a, 'XData');
y_real = get(a, 'YData');
%x_real = x_real(1:end-1);
%y_real = y_real(1:end-1);
x_real = x_real/100; % convert cm to m
y_real = y_real/100; % convert cm to m
close;
% get angle
rise = y_real(2) - y_real(1);
run = x_real(2) - x_real(1);
o = atand(rise/run); % measure launch angle (doesn't match well with slant)
% get launch velocity
disp('actual launch angle')
o
delta_t = 1/60; % based on frame rate 60
delta_d = sqrt(rise^2 + run^2); %distance traveled from initial launch
v0 = delta_d/delta_t; %cm to m

% set variables
m = .0002; % .2g to kg
d = .023; % 2.3cm diameter to m
A = pi*(d/2)^2; %area

%%%%%%% Theoretical (no drag) %%%%%%%%
dt = .001;
t(1) = 0;
vx(1) = v0*cosd(o);
vy(1) = v0*sind(o);
rx(1) = x_real(1); % Initial x from first data points
ry(1) = y_real(1); % Initial x from first data points
i = 1;
while ry(i) > 0
t(i+1) = t(i)+dt;
vx(i+1) = vx(i);
ay = -9.8; % acceleration from gravity
vy(i+1) = vy(i)+ay*dt;
rx(i+1) = rx(i)+vx(i)*dt;
ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
i = i+1;
end

rx1 = 100*rx; % save for plotting. convert to cm
ry1= 100*ry;



clear t vx vy rx ry; 
%%%%%%%% Theoretical (with air resistance) %%%%%%%%
% initial conditions
t(1) = 0;
vx(1) = v0*cosd(o);
vy(1) = v0*sind(o);
rx(1) = x_real(1); % Initial x from first data points
ry(1) = y_real(1); % Initial x from first data points
i = 1;
while ry(i) > 0
[ax, ay] = drag(vx(end),vy(end),m, d);
%ax = ...; % expression for x-component of acceleration
%ay = ...; % expression for y-component of acceleration
t(i+1) = t(i)+dt;
vx(i+1) = vx(i)+ax*dt;
vy(i+1) = vy(i)+ay*dt;
rx(i+1) = rx(i)+vx(i)*dt+ax*0.5*dt*dt;
ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
i = i+1;
end

rx2 = 100*rx; % save for plotting. convert to cm
ry2 = 100*ry;

% Plot for 1.3 platic compared to no-drag
figure
plot(rx1,ry1,rx2,ry2,100*x_real,100*y_real) % convert real to cm
legend('theoretical - no drag','theoretical - drag', 'experimental','location','southwest')
xlabel('distance (cm)')
ylabel('height (cm)')
title('plastic projectile - comparing model to experiment')

% figure
% plot(rx2,ry2,100*x_real,100*y_real) % convert real to cm
% legend('theoretical - with drag', 'experimental','location','southwest')
% xlabel('distance (cm)')
% ylabel('height (cm)')
% title('plastic projectile - comparing model to experiment')


%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%

function [ax, ay] = drag(vx_i,vy_i,m_i,d)
	g = 9.81; % gravitational acceleration (always in -y direction)
	C = 1; % Drag coefficient
	p = 1.2; % air density (kg/m^3)
	r = d/2; % diameter -> radius
	A = pi*r^2; % cross sectional area
	%o = atand(vy_i/vx_i); % angle of travel
	vi = sqrt(vx_i^2 + vy_i^2); %initial velocity
	d = 1/2*C*A*p*vi; % drag proportional to v

	ax = -d/m_i*vx_i; % negative counters direction. *v to get v^2
	ay = -g - d/m_i*vy_i; % negative counters direction. *v to get v^2
end
