% Matlab sucks
% Successorive over-relaxation

A = [2 -1 0; 
    -1 2 -1;
    0 -1 2];

c = [1;4;7];

% A is positive definite if any arbitrary vector v gives a positive value when
% put into the following form:
% vt A v (vt = v transpose)

v = [1 ;1 ;1];
vt = v.';  % alternatively transpose(v)

% for example:
vt*A*v % +2 => positive definite

% in general, if we let v = [x y z]
% matrix looks like
%                                 |x|
% |(2x - y) (-x + 2y -z) (-y +2z)||y|   = 2x^2 - xy -xy + 2y^2 -yz -yz + +2z^2
%                                 |z|
%
% = 2x^2 -2xy + y^2 -2yz 2z^2 = x^2 + (x-y)^2 + (y-z)^2 + z^2

% solve with Gauss-Seidel

[x n] = gauss_seidel(A,c);

fprintf('After %d iterations, found value of x as:', n)
x

% Use Succesive Order Relaxation to find the optimum weight

ws = [];
ns = [];
sol_real = [4.5; 8; 7.5];
for w = 0:.1:2;
	[x n] = gauss_seidel_SOR(A,c,w);
	ns = [ns n];
	ws = [ws w];
end

% plot results of SOR

plot(ws, ns)

% The optimal w appears to be around w = 1.3

%%%%%%%%%%%$$$$$$$$$$%% functions %%%%%%%%%%%%%%%%%%%%%%%%

function [x n] = gauss_seidel(A,c)
	M = tril(A);
	N = M - A;
	x = [0;0;0]; %initial guess
	change = 100;
	n = 0 ;
	while change > 1e-6
		x_old = x;
		x = inv(M)*N*x + inv(M)*c;
		change = abs(x - x_old);
		n = n + 1;
	end
	
end

function [x n] = gauss_seidel_SOR(A,c,w)
	M = tril(A);
	N = M - A;
	x = [0;0;0]; %initial guess
	change = 100;
	n = 0 ;
	while change > 1e-6
		x_old = x;
		x = inv(M)*N*x + inv(M)*c;
		x = w*x + (1-w)*x_old;
		change = abs(x - x_old);
		n = n + 1;
		if n > 500
			break;
		end
	end
	
end
