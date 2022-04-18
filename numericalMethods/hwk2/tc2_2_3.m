% tridiagonal matrix system - Thomas algorithm

A =    [1 2 0 0 0;
	-2 4 5 0 0;
	0 7 -9 2 0;
	0 0 6 3 2;
	0 0 0 5 3]

b =   [-1;
	2;
	-3;
	4;
	3]


x = thomas(A,b);
x

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
