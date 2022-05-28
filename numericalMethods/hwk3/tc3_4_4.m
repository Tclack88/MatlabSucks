% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'

clc
clear

N = 40;
x = linspace(0,20,N);
dx = x(2) - x(1)

dts = [.01 .05 .1 .2 .5 1];
figure
hold on
for dt = dts;
	O = initiate_phis(x);
	evals = get_eig_vals(O,dx)';
	X = real(evals)*dt;
	Y = imag(evals)*dt;
	plot(X,Y,'o');
end
legend_strings = "dt: " + string(dts);
legend(legend_strings,'location','southwest')
title('stablility plots for various dt time steps with RK2')


% add stability plot over this
relamdt=-3:0.1:1;
imlamdt=-2:0.1:2;
[RElamdt,IMlamdt]=meshgrid(relamdt,imlamdt);
axis square;
lamdt=RElamdt+i*IMlamdt;
sig=1+lamdt+lamdt.^2/2;
contour_levels=[1 0 0 0];
contour(RElamdt,IMlamdt,abs(sig),contour_levels,'linewidth',4);
xlabel('Re \lambda \Delta t');
ylabel('Im \lambda \Delta t');
hold on;


disp('we can see certain eigenvalues vall outside of the region of stability')
disp('Some eigenvalues for dt=1 are outside the range of stability')
disp('this makes this division unsuitable for stabile solutions')


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
	M1 = -U/(2*dx)*M1;
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
