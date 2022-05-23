% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'

N = 40;
x = linspace(0,20,N);
dx = x(2) - x(1)

dts = [.01 .05 .1 .2 .5];
figure
hold on
for dt = dts
	O = initiate_phis(x);
	evals = get_eig_vals(O,dx)';
	X = real(evals);
	Y = imag(evals);
	for x_val = X
		if -2 > x_val > 0
			disp('outside of stablility region')
			dt
			x_val
		end
	end
	plot(X,Y,'o');
end

% plot phis for when t = 0,5,10,15
% times = [0 5 10 15];
% for time = times
% 	t_i = time/dt + 1; %index with our phis
% 	plot(x,O(t_i,:))
% 	xlabel('x')
% 	ylabel('phi')
% 
% end
% The equation we are dealing with:
% dO/dt + UdO/dx = v*d2O/dx2"

function evals = get_eig_vals(O,dx)
	% set constants
	v = .025;
	U = 1;
	% Return derivatives in vector:
	O_0 = O(1);
	L = L_mat(O,dx,U,v);
	evals = eig(L);
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
