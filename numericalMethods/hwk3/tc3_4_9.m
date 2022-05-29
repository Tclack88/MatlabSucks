% Unsteady convection -- diffusion
% dO/dt + U*dO/dx = v*D2O/Dx2
% Use N=40 and dt = .1
clc
clear

dt = .001;
a = 0; % left endpoint for both t and x
b = 20;% right endpoint for both t and x
t = [a:dt:b];
N=100;
N_1 = N-1; % N-1, GLL nodes gives one more point than the number passed in
x = GLL_nodes(N_1);
x = (b-a)/2*x + (b+a)/2; % Adjust Gauss lobatto beyond the -1 +1 region
D = DerivMatrix(x,N_1);

O = zeros(length(t)-1,N);
O(1,:) = initiate_phis(x);


for n=1:length(t)-1
	k1=f(D, O(n,:));
        k2=f(D, O(n,:)+dt*k1);
	O(n+1,:) = O(n,:) + (k1/2+k2/2)*dt;
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
title('Phi at x for with dt=.001 and N=100 using spectral differentiation')
legend('time:'+string(times))
hold off


disp('graphically, with dt=.001 the results are as expected, agreeing with 4.3')

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
