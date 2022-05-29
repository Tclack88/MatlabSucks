% Unsteady convection -- diffusion
% dO/dt + U*dO/dx = v*D2O/Dx2
% Use N=100 and find stable dt values
clc
clear

N=100;
a = 0;
b = 20;
dt=.1; 
t = [a:dt:b];
N_1 = N-1; % N-1, GLL nodes gives one more point than the number passed in
x = GLL_nodes(N_1);
x = (b-a)/2*x + (b+a)/2; % Adjust Gauss lobatto beyond the -1 +1 region
D = DerivMatrix(x,N_1);

O = zeros(length(t)-1,N);
O(1,:) = initiate_phis(x);


figure
hold on
% overly plots with different dts to estimate stability
dts = [1 .1 .01 .001];
for dt = dts;
        evals = get_eig_vals(D)';
        X = real(evals)*dt;
        Y = imag(evals)*dt;
        plot(X,Y,'o');
end

% add stability plot
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
hold off;
xlabel('x')
ylabel('phi')
title('stability plot for several values of dt')
legend_strings = "dt: " + string(dts);
legend(legend_strings,'location','southwest')

disp('part of dt=.01 is within the region of stability as well as much of dt=.001')



%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%
function evals = get_eig_vals(D)
        % set constants
        v = .025;
        U = 1;
        % build matrix
	D2 = D*D;
	L = v*D2 - U*D;
        evals = eig(L);
end

function vec = f(D,O)
	% Return next set of phis
	v = .025;
	U = 1;
	O(1) = 0; % boundary condition requires this
	D2 = D*D;
	% Return derivatives in vector:
	L = v*D2 - U*D;
	vec = L*O';
	vec = vec';
end


function O = initiate_phis(x)
	% pass in array of x,
	% returns starting phis
	O = exp(-((x-3)/.5).^2); % from BC
	O(1) = 0; % from BC
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
