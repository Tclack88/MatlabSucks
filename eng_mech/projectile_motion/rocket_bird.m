clc;
clear;
m = 1.2;
m_i = 1.2;
v0 = 10;
o = 50;


%%%%%%% original no drag %%%%%%%%
dt = .001;
t(1) = 0;
vx(1) = v0*cosd(o);
vy(1) = v0*sind(o);
rx(1) = 0;
ry(1) = 20;
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
rx1 = rx;
ry1 = ry;



%%%%%%%% original with air resistance %%%%%%%%
dt = 0.001;
% initial conditions
clear t vx vy rx ry;
t(1) = 0;
vx(1) = v0*cosd(o);
vy(1) = v0*sind(o);
rx(1) = 0;
ry(1) = 20;
i = 1;
while ry(i) > 0
[ax, ay] = drag(vx(end),vy(end),m);
t(i+1) = t(i)+dt;
vx(i+1) = vx(i)+ax*dt;
vy(i+1) = vy(i)+ay*dt;
rx(i+1) = rx(i)+vx(i)*dt+ax*0.5*dt*dt;
ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
i = i+1;
end
rx2 = rx;
ry2 = ry;

clear t vx vy rx ry;

dt = 0.001;
% initial conditions
% clear t;
% clear vx;
% clear vy;
% clear rx;
% clear ry
t(1) = 0;
vx(1) = v0*cosd(o);
vy(1) = v0*sind(o);
rx(1) = 0;
ry(1) = 20;
i = 1;
while ry(i) > 0
[ax, ay] = drag(vx(end),vy(end),m);
t(i+1) = t(i)+dt;
if mod(t(i),.1) == 0 || mod(t(i),.1) <= 1.0000e-03
	% matlab, a terrible "programming language", didn't see ANY of the .1
	% moded values after the first, so the extra "or" condition is necessary
	% for it to cause mass_ejection to occur
	[vx_f,vy_f,m_f] = mass_ejection(vx(i),vy(i),50,m,.01*m_i);
	vx(i+1) = vx_f;
	vy(i+1) = vy_f;
	m = m_f; % set new mass
else
	vx(i+1) = vx(i)+ax*dt;
	vy(i+1) = vy(i)+ay*dt;
end
rx(i+1) = rx(i)+vx(i)*dt+ax*0.5*dt*dt;
ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
i = i+1;
end
rx3 = rx;
ry3 = ry;
figure
plot(rx1, ry1, rx3,ry3) % original (no drag) + rocket bird)
title('Trajectory comparing original to mass ejection and drag')
xlabel('x (cm)')
ylabel('y (cm)')
legend('original', 'rocket bird')

%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%

function [ax, ay] = drag(vx_i,vy_i,m)
	g = 9.81; % gravitational acceleration (always in -y direction)
	C = 1; % Drag coefficient
	p = 1.2; % air density (kg/m^3)
	r = .15; % sphere radius (m)
	A = pi*r^2; % cross sectional area

	o = atand(vy_i/vx_i); % angle of travel
	vi = sqrt(vx_i^2 + vy_i^2); %initial velocity
	d = 1/2*C*A*p*vi; % drag proportional to v

	ax = -d/m*vx_i; % negative counters direction. *v to get v^2
	ay = -g - d/m*vy_i; % negative counters direction. *v to get v^2
end


function [vx_f,vy_f,m_f] = mass_ejection(vx_i,vy_i,v_eb,m_i,m_e)
	% Inputs:
	% vx_i: Velocity component in x direction before ejection
	% vy_i: Velocity component in y direction before ejection
	% v_eb: Ejection velocity of the mass relative to the bird
	% m_i : Mass of the bird before ejection
	% m_e : Ejection mass
	% Outputs:
	% vx_f: Velocity component in x direction after ejection
	% vy_f: Velocity component in y direction after ejection
	% m_f : Mass of the bird after ejection
	m_f = m_i - m_e;
	vi = sqrt(vx_i^2 + vy_i^2); % instantaneous magnitude of velocity
	o = atand(vy_i/vx_i); % instantaneous angle of travel (radians)
	ve = vi - v_eb;
	vf = (m_i*vi - m_e*ve)/(m_i - m_e);
	vx_f = vf*cosd(o);
	vy_f = vf*sind(o);
end
