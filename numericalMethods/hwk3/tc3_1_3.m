% Thomas algorithm to compute Boundary value ODE of y'' + y' - 6y = 15sin(12x)

deltas = [.2 .1 .01];
left = -1; right =2;
hold on
for delta = deltas;
	xp = [left:delta:right];
	D0 = YMatrix(xp);
	D1 = D1Matrix(xp);
	D2 = D2Matrix(xp);
	b = 15*sin(12*xp);
	%b(1) = left;
	%b(end) = right;
	b = b'; %transpose for correct shape
	D = D2 + D1 - D0;
	size(b);
	size(D);
	sol = linsolve(D,b);
	plot(xp,sol)
end
hold off
title('ODE plot of y''''+y''-6y = 15sin(12x) at different intervals');
legendStrings = "delta:" + string(deltas);
legend(legendStrings,'location','northwest')

function D = D2Matrix(xp)
	N = length(xp);
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


function D = D1Matrix(xp)
	N = length(xp);
        D = zeros(N);
        D(:,1) = 1;
        D(:,2) = 0;
        D(:,3) = 1;
        for i = 3:N-1;
                D(i,:) = circshift(D(i,:),i-2);
        end
        D(end,:) = circshift(D(end,:),N-3);
        delta = xp(2) - xp(1); % same as interval linspace(a,b,N+1)
        D = D/(2*delta);
end
	
function D = YMatrix(xp)
	N = length(xp);
	D = 6*eye(N);
end
