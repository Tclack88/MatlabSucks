delta = 1.3;
l = 130;
E = 3.1e7;
S = 1100;
I = 630;
xs = [0:delta:l];
N = length(xs);
qs = [10 15 25];

y_lim = 1/285; % Establish upper bound to check against
% By symmetry, max deflection occurs right at center (x=65in), we can slowely
% increment q until we reach the limit y_max set above
q = 0;
y_max = 0;
while y_max < y_lim
	q = q+1;
	Q = q/(2*E*I)*xs.*(xs-l);
	Q = Q'; % transpose fit matrix
	A = D2Matrix(xs) - S/(E*I)*eye(length(xs));
	A(1) = 1;
	A(end) = 1;
	sol = linsolve(A,Q);
	y_max = max(sol);
end
q_max = q-1; % subtract 1 because last q was above our limit
sprintf('maximum allowed loading (to nearest lb) found to be %d lbs/in',q_max)

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

