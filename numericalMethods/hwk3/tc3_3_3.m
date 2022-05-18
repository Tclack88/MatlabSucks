% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'
% Compare with different values for A


dt = .01;
t = [0:dt:50];

Os = [0 10 30 50 90];
for o = Os
	x = zeros(length(t),2);
	x(1,:) = [0,0];
	for n=1:length(t)-1
		k1=f(t(n), x(n,:),o)';
		k2=f(t(n)+dt,x(n,:)+dt*k1,o)';
		x(n+1,:)=x(n,:) + (k1/2+k2/2)*dt;
	end
	plot(t,x(:,2))
	hold on
end

title("RK2 for system of ODEs x'' = 50Vcos(o) + 6Asin(x/10) - 60x'")
subtitle("comparing different values of o (angle wind blows)")
legend_string = "angle= " + string(Os);
legend(legend_string,'location','southeast')

disp('As o (wind angle) increases from 0 to 90, the velocity of the boat drops. Ths makes sense as the angle is in the direction of the boat')
disp('a smaller angle means more wind is "boosting" the boat along. When o = 90, the wind is blowing directly against')
disp('the boat and with these particular parameters, the boat does not move (vx = 0)')
% The equation can be rearranged as:
% x'' = 1/2Vcos(o) + 3/50Asin(x/10) - 3/5x'
% Here, we can reduce to a system of equations:
% x1 = x      x2 = x'
% so we have:
% x1' = x2
% x2' = 1/2Vcos(o) + 3/50Asin(x1/10) - 3/5(x2)

function vec = f(t,x,o)
	% set constants
	V = 14;
	A = 6;
	% Return derivatices in vector:
	% [x1' x2']
	vec = [x(2); (1/2*V*cosd(o) + 3/50*A*sin(x(1)/10)) - 3/5*x(2) ];
end
