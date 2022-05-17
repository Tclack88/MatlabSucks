% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'


dt = .1;
t = [0:dt:50];
x = zeros(length(t),2);
x(1,:) = [0,0];

for n=1:length(t)-1
	k1=f(t(n), x(n,:))';
	k2=f(t(n)+dt,x(n,:)+dt*k1)';
	x(n+1,:)=x(n,:) + (k1/2+k2/2)*dt;
end;

[tmat,xmat]=ode23(@f,[0 50];[0 0])

plot(t,x(:,2))
hold on
plot(tmat,xmat(:,2))

% The equation can be rearranged as:
% x'' = 1/2Vcos(o) + 3/50Asin(x/100) - 3/5x'
% Here, we can reduce to a system of equations:
% u = x  u'=v
% v = x' v'= 1/2Vcos(o) + 3/50Asin(u/100) - 3/5v
%   0          1
% -3/5   1/2Vcos(o) + 3/50Asin(u/100) 

function vec = f(t,x)
	% set constants
	V = 14;
	o = 30;
	A = 3;
	% system of equations
	vec = [x(2); (1/2*V*cosd(o) + 3/50*A*sin(x(1)/100))];
end
