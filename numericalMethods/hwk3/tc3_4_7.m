% Unsteady convection -- diffusion
% dO/dt + U*dO/dx = v*D2O/Dx2
% Use N=40 and dt = .1
clc
clear

dt = .1;
a = 0; % left endpoint for both t and x
b = 20;% right endpoint for both t and x
t = [a:dt:b];
N=40;
N_1 = N-1; % N-1, GLL nodes gives one more point than the number passed in
%x = linspace(0,20,N);
x = GLL_nodes(N_1);
%[w1, x1] = gauss_legendre_points(N);
%[x,p] = lglnodes(N_1);
%x = x(end:-1:1); % reverse if using lglnodes because it gives it backwards
x = (b-a)/2*x + (b+a)/2; % Adjust Gauss lobatto beyond the -1 +1 region
D = DerivMatrix(x,N_1);
D

O = zeros(length(t)-1,N);
O(1,:) = initiate_phis(x);


for n=1:length(t)-1
	k1=f(D, O(n,:));
        k2=f(D, O(n,:)+dt*k1);
	O(n+1,:) = O(n,:) + (k1/2+k2/2)*dt;
	%O(n+1,:) = f(D, O(n,:))*dt;
end;

% plot phis for when t = 0,5,10,15
times = [0 5 10 15];
figure
hold on
for time = times
	t_i = time/dt + 1; %index with our phis
	plot(x,O(t_i,:))
end
xlabel('x')
ylabel('phi')
%title('Phi at x for with dt=.1 and N=40 using spectral differentiation')
legend('time:'+string(times))
hold off

O

function vec = f(D,O)
	% Return next set of phis
	v = .025;
	U = 1;
	O(1) = 0; % boundary condition requires this
	D2 = D*D;
	% Return derivatives in vector:
	L = .01*v*D2 - .1*U*D/2; % scale D by (b-a)/2 and D2 by the square
	%L = v*D2 - U*D/2;
	vec = L*O';
	vec = vec';
end


function O = initiate_phis(x)
	% pass in array of x,
	% returns starting phis
	O = exp(-((x-3)/.5).^2); % from BC
	O(1) = 0; % from BC
end

function [w,x] = gauss_legendre_points(k)
    syms t
    x = double(vpasolve(legendreP(k,t) == 0));
    w = 2*(1-x.^2)./((k+1)^2.*(legendreP(k+1,x).^2));
end

function xp = GLL_nodes(N)
        % provided in solutions workshop 8
        syms x
        Lo = (1-x^2)*diff(legendreP(N,x),x);
        xp = double(vpasolve(Lo == 0));
end

function D = DerivMatrix(x,n)
        % Spectral differentiation
        % Find D matrix (given in lecture)
        D = zeros(n+1,n+1);
        for i = 1:n+1;
                num(i)=1.0;
                for k=1:n+1
                        if k~=i
                                num(i) = num(i)*(x(i)-x(k));
                        end
                end
        end
        for j=1:n+1;
                den(j)=1.0;
                for k =1:n+1
                        if k~=j
                                den(j) = den(j)*(x(j)-x(k));
                        end
                end
        end
        for i=1:n+1
                for j=1:n+1
                        if i~=j
                                D(i,j)=num(i)/(den(j)*(x(i)-x(j)));
                        end
                end
        end
        for i=1:n+1
                for k=1:n+1
                        if i~=k
                                D(i,i) = D(i,i) + 1./(x(i)-x(k));
                        end
                end
        end
end

function [x,P]=lglnodes(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% lglnodes.m
%
% Computes the Legendre-Gauss-Lobatto nodes, weights and the LGL Vandermonde
% matrix. The LGL nodes are the zeros of (1-x^2)*P'_N(x). Useful for numerical
% integration and spectral methods.
%
% Reference on LGL nodes and weights:
%   C. Canuto, M. Y. Hussaini, A. Quarteroni, T. A. Tang, "Spectral Methods
%   in Fluid Dynamics," Section 2.3. Springer-Verlag 1987
%
% Written by Greg von Winckel - 04/17/2004
% Contact: gregvw@chtm.unm.edu
% This MATLAB Function is taken from MATLAB Central File Exchange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncation + 1
N1=N+1;
% Use the Chebyshev-Gauss-Lobatto nodes as the first guess
x=cos(pi*(0:N)/N)';
% The Legendre Vandermonde Matrix
P=zeros(N1,N1);
% Compute P_(N) using the recursion relation
% Compute its first and second derivatives and
% update x using the Newton-Raphson method.
xold=2;
while max(abs(x-xold))>eps
    xold=x;

    P(:,1)=1;    P(:,2)=x;

    for k=2:N
        P(:,k+1)=( (2*k-1)*x.*P(:,k)-(k-1)*P(:,k-1) )/k;
    end
end
end
