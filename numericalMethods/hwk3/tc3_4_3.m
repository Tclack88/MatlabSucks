% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'

% N = 6;
% x = linspace(0,20,N);
% L(x,1,1)
% return;

dt = .1;
t = [0:dt:20];
N = 20;
x = linspace(0,20,N);
dx = x(2) - x(1)

O = zeros(length(t)-1,N);
O(1,:) = initiate_phis(x);


for n=1:length(t)-1
	k1=f(t(n), O(n,:),dx)';
	k2=f(t(n)+dt,O(n,:)+dt*k1,dx)';
	O(n+1,:)=O(n,:) + (k1/2+k2/2)*dt;
end;

% plot phis for when t = 0,5,10,15
times = [0 5 10 15];
for time = times
	t_i = time/dt + 1; %index with our phis
	plot(x,O(t_i,:))
	xlabel('x')
	ylabel('phi')

end
% Compare RK2 with a built in
% [tmat,xmat]=ode23(@f,[0 20],O(1,:));
% hold on
% plot(tmat,xmat(:,2))
% title("RK2 for system of ODEs dO/dt + UdO/dx = v*d2O/dx2")
% legend('RK2 approximation','ode23 - Matlab builtin solver','location','southeast')
% 
% disp('the matlab built-in matches very well with the RK2 solution, so in a handy-wavy way, the chosen delta t of .1 is a good choice')

% The equation we are dealing with::
% dO/dt + UdO/dx = v*d2O/dx2"

function vec = f(t,O,dx)
	% set constants
	v = .025;
	U = 1;
	% Return derivatives in vector:
	O_0 = O(1);
	L = L_mat(O,dx,U,v);
	P = zeros(length(O)-1,1)';
	P(1) = U/(2*dx) + v/dx^2;
	vec = L*O(2:end)' + P'*O_0; % shouldn't his be 5x1?
	vec = [O_0; vec];
end

function L = L_mat(O,dx,U,v)
	N = length(O)-1; %eg:6 Phi vals, 1st and last go into P. need 5x5 matrix
	% Make M1
	M1 = zeros(N);
	M1(2:end-1,1) = -1;
        M1(2:end-1,3) = 1;
        for i = 3:N-1;
        	M1(i,:) = circshift(M1(i,:),i-2);
	end
	M1(1,2) = 1;
	M1(end, end-1) = -2;
	M1(end,end) = 2;
	M1 = -U/(2*dx);
	% Make M2
	M2 = zeros(N);
	M2(2:end,1) = 1;
        M2(2:end,2) = -2;
        M2(2:end,3) = 1;
        for i = 3:N-1;
        	M2(i,:) = circshift(M2(i,:),i-2);
	end
	M2(end,:) = circshift(M2(end,:),N-3);
	M2(1,1) = -2;
	M2(1,2) = 1;
	M2 = v/dx^2*M2;
	L = M1 + M2;
end	


function O = initiate_phis(x)
	% pass in array of x,
	% returns starting phis
	O = exp(-((x-3)/.5).^2); % from BC
	O(1) = 0; % from BC
end
