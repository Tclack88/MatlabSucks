% Partial double derivatives of a 2D function

% Create the D2 Matrix (matrix that gives 2nd derivatives)

a = -1; % left boundary
b = 1; % right boundary
N = 50; % number of points/size of matrix
xp = linspace(a,b,N)'; % need a column matrix
u = xp.^7;

D2 = SecondDerivMatrix(N,xp);
D2u = D2*u; % numerical approximation solution
u_r = 42*(xp.^5); % exact solution


% plot of 2nd derivative of x^7
plot(xp,D2u,xp,u_r,'.m');
title('Plot of derivative of x^7')
xlabel('x')
ylabel('x^7')
legend('numerical solution', 'exact solution')

disp('50 points shown, but after about 15-20 points, the numerical solution lines up nicely with the approximate solution')



% Now we check partial 2nd derivatives to plot a 2d function (in 3D)
N = 50; %number of points/size of matrix
x = linspace(-1,1,N);
y = linspace(0,2,N);

[X, Y] = meshgrid(x,y);
u = (sin(Y).*(X.^5))'; % function points (transposed to column for matrix mult)

uxx = 20*sin(Y).*(X.^3);
uyy = -sin(Y).*(X.^5);

figure
surf(uxx + uyy, X,Y);
title('Actual solution plot of sin(y)x^5 2nd partials')



% Numerical approximation of the 2nd partial derivative
D = SecondDerivMatrix(N,x);
I = eye(N);

D2x = kron(I,D); % kronecker tensor product
D2y = kron(D,I);

u = u';
u = u(:);

uxx = D2x*u;
uyy = D2y*u;

% reshape for plotting
% My first instinct was to transpose uxx and uyy but that was wrong
uxx = reshape(uxx,N,N);
uyy = reshape(uyy,N,N);
figure
surf(uxx + uyy,X,Y);
title('Numerical methods approach plot of sin(y)x^5 2nd partials')



%%%%%% functions %%%%%%%%
function D = SecondDerivMatrix(N,xp)
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

