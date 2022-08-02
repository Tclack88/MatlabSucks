m=25;
N = 6


x = linspace(0,1);
c = zeros(1,N);

O = zeros(size(x));
O(1) = 0;
O(end) = 1;



D1 = D1Matrix(N,x)
D2 = D2Matrix(x)


A = D2 + (2/x)*D1 - (m/x)*eye(N); % matrix representic system
% this doesn't work because x is a factor out front, 
%  :[

sol = linsolve(A,x)
A = D2Matrix(xs) - S/(E*I)*eye(length(xs));
sol = linsolve(A,c);


function D = D1Matrix(N,xs)
        % specific parameters

        delta = xs(2) - xs(1); % same as interval linspace(a,b,N+1)

        D = zeros(N);
        D(2:end-1,1) = -1;
        D(2:end-1,2) = 0;
        D(2:end-1,3) = 1;
        for i = 3:N-1;
                D(i,:) = circshift(D(i,:),i-2);
        end
	D = D/(2*delta)
        %D(1) = 1;
        %D(end) = 1;

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


