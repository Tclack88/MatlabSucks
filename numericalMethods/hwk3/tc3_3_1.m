% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'


dt = .01;
t = [0:dt:50];
x = zeros(length(t),2);
x(1,:) = [0,0];

for n=1:length(t)-1
	k1=f(t(n), x(n,:))';
	k2=f(t(n)+dt,x(n,:)+dt*k1)';
	x(n+1,:)=x(n,:) + (k1/2+k2/2)*dt;
end;

plot(t,x(:,2))
% Compare RK2 with a built in
[tmat,xmat]=ode23(@f,[0 50],[0 0]);
hold on
plot(tmat,xmat(:,2))
title("RK2 for system of ODEs x'' = 50Vcos(o) + 6Asin(x/10) - 60x'")
legend('RK2 approximation','ode23 - Matlab builtin solver','location','southeast')

disp('the matlab built-in matches very well with the RK2 solution, so in a handy-wavy way, the chosen delta t of .1 is a good choice')

% The equation can be rearranged as:
% x'' = 1/2Vcos(o) + 3/50Asin(x/10) - 3/5x'
% Here, we can reduce to a system of equations:
% x1 = x      x2 = x'
% so we have:
% x1' = x2
% x2' = 1/2Vcos(o) + 3/50Asin(x1/10) - 3/5(x2)

function vec = f(t,x)
	% set constants
	V = 14;
	o = 30;
	A = 3;
	% Return derivatices in vector:
	% [x1' x2']
	vec = [x(2); (1/2*V*cosd(o) + 3/50*A*sin(x(1)/10)) - 3/5*x(2) ];
end
