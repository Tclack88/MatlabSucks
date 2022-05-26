% evolution of N vortices in an inciscid fluid

dt = .01;
t = [0:dt:200];
size(t);
200/dt;
x_space = linspace(-5,5);
y_space = linspace(-5,5);
v1 = [1; 0]; % random starting point
v2 = [-1; 0]; % random starting point
N = 2;

vs = zeros(4,length(t));
vs(:,1) = [v1;v2];
for n = 1:length(t)-1;
	dvdt = set_2_vortex_eqs(vs(:,n),N);
	vs(:,n+1) = vs(:,n) + dvdt*dt;
end

hold on
plot(vs(1,:),vs(2,:))
plot(vs(3,:),vs(4,:))
legend('vortex1','vortex2')
xlabel('x position')
ylabel('y position')


function dvdt = set_2_vortex_eqs(v,N);
	%w = 1; % ignore for simplicity, set to 1
	M = 	[0 -1 0 1; %multiply by 1 to keep same factor 1/2pi
	         1 0 -1 0;
	 	 0 1 0 -1; %multiply by 1 to keep same factor 1/2pi
		 -1 0 1 0];
	r = sqrt((v(1)-v(3))^2 + (v(2)-v(4))^2);
	M = 1/(2*pi)*M/(r^2);
	dvdt = M*v;
end
