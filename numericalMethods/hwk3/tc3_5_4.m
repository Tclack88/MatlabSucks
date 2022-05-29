% Poisson solver for gauss-lobatto nodes, and exact solution.
clear
clc

n = 20;

method = 'Even';
[X,Y,uzz] =  Poisson_solve(@f,n,method) 
figure
surf(X,Y, uzz);
title({'Poisson solve N=',string(n),method})


figure
method = 'GLo';
[X,Y,uzz] =  Poisson_solve(@f,n,method) 
surf(X,Y, uzz);
title({'Poisson solve N=',string(n),method})

%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%

function u = f(X,Y)
	u = 10*sin(8*X.*(Y-1));
end

function [X,Y,uzz] =  Poisson_solve(f,N,method) 
	if strcmp(method,'Even');
		x = linspace(-1,1,N);
		y = linspace(-1,1,N);
		D2 = SecondDerivMatrix(x,N);
	else strcmp(method,'GLo');
		x = GLL_nodes(N-1);
		y = GLL_nodes(N-1);
		D = DerivMatrix(x,N-1);
		D2 = D*D;
	end
	[X,Y] = meshgrid(x,y);
	u = f(X,Y)
	% Set boundary conditions
	u(1,:) = 0;
	u(end,:) = 0;
	u(:,1) = 0;
	u(:,end) = 0;
	I = eye(N);
	D2x = kron(I,D2);
	D2y = kron(D2,I);
	u = u(:);
	uzz = (D2x + D2y)*u;
	uzz = reshape(uzz,N,N);
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

function D = SecondDerivMatrix(xp,N)
        D = zeros(N);
        D(:,1) = 1;
        D(:,2) = -2;
        D(:,3) = 1;
        for i = 3:N-1;
                D(i,:) = circshift(D(i,:),i-2);
        end
        D(end,:) = circshift(D(end,:),N-3);
        delta = xp(2) - xp(1); % same as interval linspace(a,b,N+1)
        D = D/(delta^2);
end

function xp = GLL_nodes(N)
	% provided in solutions workshop 8
	syms x
	Lo = (1-x^2)*diff(legendreP(N,x),x);
	xp = double(vpasolve(Lo == 0));
end
