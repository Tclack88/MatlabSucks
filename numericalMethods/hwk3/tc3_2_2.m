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
	A = SetMatrix(N,xs);
	sol = thomas(A,Q);
	plot(xs,sol);
end
hold off
title('Beam deflection ODE solution (y'''' = S/(EI) y + q/(2EI)x(x-l))')
legendStrings = "q = " + string(qs);
legend(legendStrings)

format short
deltas = [1 .1 .01];
% Compare the mean squared difference between the max value of the original
% (1.3 inch) gradations to some smaller choices
for delta = deltas
	xs = [0:delta:l];
	N = length(xs);
	Q = q/(2*E*I)*xs.*(xs-l);
	Q = Q'; % transpose fit matrix
	A = SetMatrix(N,xs);
	sol_new = thomas(A,Q);
	err = rms(max(sol),max(sol_new));
	sprintf('using delta = %.2f in. to 1.3 in. gradations is: %.10f %% difference, so 1.3 is a very good choice', delta,err)
end

function err = rms(a1, a2)
        % estimate error with root mean square between 2 arrays
        err = sqrt(sum((abs(a1 - a2)).^2)/length(a1));
end


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


