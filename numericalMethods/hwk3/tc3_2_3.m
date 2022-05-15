delta = 1.3;
l = 130;
E = 3.1e7;
S = 1100;
I = 630;
xs = [0:delta:l];
N = length(xs);
qs = [10 15 25];

hold on
for q = qs;
	Q = q/(2*E*I)*xs.*(xs-l);
	Q = Q'; % transpose fit matrix
	A = D2Matrix(xs) - S/(E*I)*eye(length(xs));
	A(1) = 1;
	A(end) = 1;
	sol = linsolve(A,Q);
	plot(xs,sol);
end
hold off
title('Beam deflection ODE solution (y'''' = S/(EI) y + q/(2EI)x(x-l))')
legendStrings = "q = " + string(qs);
legend(legendStrings)

function D = SetMatrix(N,xs)
	% specific parameters
	E = 3.1e7;
	S = 1100;
	I = 630;

        delta = xs(2) - xs(1); % same as interval linspace(a,b,N+1)
        a = 1/delta^2;
        B = -2/delta^2 - S/(E*I);
        y = 1/delta^2;

        D = zeros(N);
        D(2:end-1,1) = a;
        D(2:end-1,2) = B;
        D(2:end-1,3) = y;
        for i = 3:N-1;
                D(i,:) = circshift(D(i,:),i-2);
        end
        D(1) = 1;
        D(end) = 1;

end

function D = D2Matrix(xp)
        N = length(xp);
        D = zeros(N);
        D(2:end-1,1) = 1;
        D(2:end-1,2) = -2;
        D(2:end-1,3) = 1;
        for i = 3:N-1;
                D(i,:) = circshift(D(i,:),i-2);
        end
        delta = xp(2) - xp(1); % same as interval linspace(a,b,N+1)
        D = D/(delta^2);
end

