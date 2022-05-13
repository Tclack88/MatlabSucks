% Thomas algorithm to compute Boundary value ODE of y'' + y' - 6y = 15sin(12x)

deltas = [.2 .1 .01];
left = -1; right =2;
hold on
for delta = deltas;
	xp = [left:delta:right];
	N = length(xp);
	D = SetMatrix(N,xp);
	b = 15*sin(12*xp);
	b(1) = left;
	b(end) = right;
	sol = thomas(D,b);
	plot(xp,sol)
end
hold off
title('ODE plot of y''''+y''-6y = 15sin(12x) at different intervals');
legendStrings = "delta:" + string(deltas);
legend(legendStrings,'location','northwest')


function D = SetMatrix(N,xp)
        delta = xp(2) - xp(1); % same as interval linspace(a,b,N+1)
	a = 1/delta^2 - 1/(2*delta);
	B = -2/delta^2 - 6;
	y = 1/delta^2 + 1/(2*delta);

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



function x =  thomas(A,b)
        % parse Matrix for elements of thomas algorithm
        B = diag(A);
        Q = b';
        a = zeros(size(B)-1);
        y = zeros(size(B)-1);
        for i = 1:size(a);
                a(i) = A(i+1,i);
                y(i) = A(i,i+1);
        end
        a = [0 a]; % there is no alpha 1, set to zero
        y = [y 0]; % there is no gamma n, set to zero
        % set *star values (must move forward in loop as ith element depends on
        % some values from i-1th element)
        B_s = zeros(size(B));
        Q_s = zeros(size(Q));
        B_s(1) = B(1);
        Q_s(1) = Q(1);
        for i = 2:length(B)
                B_s(i) = B(i) - a(i)*y(i-1)/B_s(i-1);
                Q_s(i) = Q(i) - a(i)*Q_s(i-1)/B_s(i-1);
        end
        % Determine solutions (move in reverse loop as ith element depends on
        % some values from i+1th element)
        x = zeros(size(Q));
        x(end) = Q_s(end)/B_s(end);
        for i = length(x)-1:-1:1; % set middle values
                x(i) = (Q_s(i) - y(i)*x(i+1))/B_s(i);
        end
end

