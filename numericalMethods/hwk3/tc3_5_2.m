%example from workshop 8.7

n = 40;
x = linspace(-1,1,n)';
y = linspace(-1,1,n)';
% x = GLL_nodes(n-1);
% y = GLL_nodes(n-1);
[X,Y] = meshgrid(x,y);

u = sin(pi*X).*cos(pi*Y);
size(u)
% numerical approximation to 2nd partial derivative
D2 = SecondDerivMatrix(x,n);
D2u = D2*u;
I = eye(n);

D2x = kron(I,D2);
D2y = kron(D2,I);

u=u';
u = u(:);

uxx = D2x*u;
uyy = D2y*u;

uxx = reshape(uxx,n,n);
uyy = reshape(uyy,n,n);
figure
surf(uxx + uyy,X,Y);

% compare to exact solution
figure
u_r = -2*pi^2*sin(pi*X).*cos(pi*Y);
surf(u_r,X,Y);




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
