% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'
% Compare with different values for A


dt = .01;
t = [0:dt:50];

As = [.1 1 6];
for A = As
	x = zeros(length(t),2);
	x(1,:) = [0,0];
	for n=1:length(t)-1
		k1=f(t(n), x(n,:),A)';
		k2=f(t(n)+dt,x(n,:)+dt*k1,A)';
		x(n+1,:)=x(n,:) + (k1/2+k2/2)*dt;
	end
	plot(t,x(:,2))
	hold on
end

title("RK2 for system of ODEs x'' = 50Vcos(o) + 6Asin(x/10) - 60x'")
subtitle("comparing different values of A (wave amplitude)")
legend_string = "A= " + string(As);
legend(legend_string,'location','southeast')

disp('As A (amplitude of the waves) increases, the velocity of the boat becomes more varied. This makes sense as the waves ')
disp('will add more resistance or add height, increasing the path length, which deters the horizontal velocity of the boat')
% The equation can be rearranged as:
% x'' = 1/2Vcos(o) + 3/50Asin(x/10) - 3/5x'
% Here, we can reduce to a system of equations:
% x1 = x      x2 = x'
% so we have:
% x1' = x2
% x2' = 1/2Vcos(o) + 3/50Asin(x1/10) - 3/5(x2)

function vec = f(t,x,A)
	% set constants
	V = 14;
	o = 30;
	% Return derivatices in vector:
	% [x1' x2']
	vec = [x(2); (1/2*V*cosd(o) + 3/50*A*sin(x(1)/10)) - 3/5*x(2) ];
end
