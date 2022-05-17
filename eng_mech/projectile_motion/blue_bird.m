clc;
clear;
m = 1.2;
v0 = 10;
o = 45;
dt = .001;

%%%%%%% initial launch %%%%%%%%
[rx, ry, vx, vy, t] = projectile(dt,v0,o);

% t(1) = 0;
% vx(1) = v0*cosd(o);
% vy(1) = v0*sind(o);
% rx(1) = 0;
% ry(1) = 20;
% i = 1;
% while t(i) < .5
% t(i+1) = t(i)+dt;
% vx(i+1) = vx(i);
% ay = -9.8; % acceleration from gravity
% vy(i+1) = vy(i)+ay*dt;
% rx(i+1) = rx(i)+vx(i)*dt;
% ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
% i = i+1;
% end

% method 1 (find in n-t then rotate)
% o = atand(vy(end)/vx(end));
% vx0 = sqrt(vy(end)^2 +vx(end)^2);
% vy0 = 0;
% vx2 = vx0;
% vy2 = 0;
% A = [1 1;tand(30) tand(-30)];
% b = [3*vx0 - vx2;3*vy0 - vy2];
% disp('answer:')
% vx2
% x = linsolve(A,b);
% vx1_0 = x(1); 
% vy1_0 = tand(30)*vx1_0;
% vx3_0 = x(2);
% vy3_0 = tand(-30)*vx3_0;
% R = [cosd(o) -sind(o);sind(o) cosd(o)]; %rotation matrix
% vx1_0
% vy1_0
% v1 = R*[vx1_0;vy1_0]
% v2 = R*[vx2;vy2]
% vx3_0
% vy3_0
% v3 = R*[vx3_0;vy3_0]



% method 2
% shorter knowing they are the same velocity in x
o = atand(vy(end)/vx(end));
R = [cosd(o) -sind(o);sind(o) cosd(o)]; %rotation matrix
vx0 = sqrt(vy(end)^2 +vx(end)^2);
vy1_0 = tand(30)*vx0;
vy3_0 = tand(-30)*vx0;
v1 = norm(R*[vx0;vy1_0]);
v2 = norm(R*[vx0;0]);
v3 = norm(R*[vx0;vy3_0]);
disp('method 2 answers:')
% v1
% v2
% v3
% setup for plotting other paths
%[rx1, ry1, vx1, vy1, t1] = projectile(t(end),dt,v1,30);
[rx1, ry1, vx1, vy1, t1] = projectile2(t(end),rx(end),ry(end),dt,v0,o+30);
[rx2, ry2, vx2, vy2, t2] = projectile2(t(end),rx(end),ry(end),dt,v0,o);
[rx3, ry3, vx3, vy3, t3] = projectile2(t(end),rx(end),ry(end),dt,v0,o-30);
hold on
plot(rx, ry)
plot(rx1, ry1)
plot(rx2, ry2)
plot(rx3, ry3)
xlabel('x (cm)')
ylabel('y (cm)')
hold off

% Energy for split
KE0 = KE(m,vx0);
KE1 = KE(m/3,v1);
KE2 = KE(m/3,v2);
KE3 = KE(m/3,v3);
KEf = KE1 + KE2 + KE3;
delta_KE = KEf - KE0;
disp('Energy required for split:')
delta_KE

function [rx, ry, vx, vy, t] = projectile(dt,v0,o)
	t(1) = 0;
	vx(1) = v0*cosd(o);
	vy(1) = v0*sind(o);
	rx(1) = 0;
	ry(1) = 20;
	i = 1;
	while t(i) < .5
		t(i+1) = t(i)+dt;
		vx(i+1) = vx(i);
		ay = -9.8; % acceleration from gravity
		vy(i+1) = vy(i)+ay*dt;
		rx(i+1) = rx(i)+vx(i)*dt;
		ry(i+1) = ry(i)+vy(i)*dt+ay*0.5*dt*dt;
		i = i+1;
	end
end

function [rx, ry, vx, vy, t] = projectile2(t0,x0,y0,dt,v0,o)
	t(1) = t0;
	vx(1) = v0*cosd(o);
	vy(1) = v0*sind(o);
	rx(1) = x0;
	ry(1) = y0;
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
end

function E = KE(m,v)
	E = 1/2*m*v^2;
end
