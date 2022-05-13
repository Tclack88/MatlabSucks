m = 1.2;
v0 = 25;
o = 45;


%%%%%%% no drag %%%%%%%%
dt = .001;
t(1) = 0;
vx(1) = v0*cos(o);
vy(1) = v0*sin(o);
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
plot(t,ry)
title('no air resistance')

clear;
clc;
m = 1.2;
v0 = 25;
o = 45;
%%%%%%%% with air resistance %%%%%%%%
dt = 0.001;
% initial conditions
t(1) = 0;
vx(1) = v0*cos(o);
vy(1) = v0*sin(o);
rx(1) = 0;
ry(1) = 20;
i = 1;
while ry(i) > 0
[ax, ay] = drag(vx(end),vy(end),m);
%ax = ...; % expression for x-component of acceleration
%ay = ...; % expression for y-component of acceleration
t(i+1) = t(i)+dt;
vx(i+1) = vx(i)+ax*dt;
vy(i+1) = vy(i)+ay*dt;
rx(i+1) = rx(i)+vx(i)*dt+ax*0.5*dt*dt;
ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
i = i+1;
end
figure
plot(t,ry)
title('with air resistance')

clear;
clc;
m = 1.2;
m_i = 1.2;
v0 = 25;
o = 45;
%%%%%%%% with air resistance %%%%%%%%
dt = 0.001;
% initial conditions
t(1) = 0;
vx(1) = v0*cos(o);
vy(1) = v0*sin(o);
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
	disp('new mass')
	m_f
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
figure
plot(t,ry)
title('with mass ejection')
disp('done')


%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%

function [ax, ay] = drag(vx_i,vy_i,m_i)
	g = 9.81; % gravitational acceleration (always in -y direction)
	C = 1; % Drag coefficient
	p = 1.2; % air density (kg/m^3)
	r = .15; % sphere radius (m)
	A = pi*r^2; % cross sectional area

	o = atan(vy_i/vx_i); % angle of travel
	vi = sqrt(vx_i^2 + vy_i^2); %initial velocity
	d = 1/2*C*A*p*vi; % drag proportional to v

	ax = -d/m_i*vx_i; % negative counters direction. *v to get v^2
	ay = -g - d/m_i*vy_i; % negative counters direction. *v to get v^2
end


function [vx_f,vy_f,m_f] = mass_ejection(vx_i,vy_i,v_eb,m_i,m_e)
	disp('mass ejection happening')
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
	o = atan(vy_i/vx_i); % instantaneous angle of travel (radians)
	vb = v_eb*m_e/(m_f); % |v| now of bird (previously 0) in bird frame
	% get x and y velocities now in bird ref. frame
	vx_b = vb*cos(o);
	vy_b = vb*sin(o);
	% get these veloxities in fixed frame
	vx_f = vx_i + vx_b;
	vy_f = vy_i + vy_b;
	
end
