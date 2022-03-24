% Gaussian Elimination. Works for an NxN matrix (augmented as well)
A0 = [1 2 3 4 8;
	2 4 8 9 12; 
	8 22 4 2 17;
	5 26 17 28 19;
	25 2 8 18 4];

b0 = [3;4;8;10;20];
b = [3;4;8;10;20];

A = [A0 b] 

x = [1;2;3;4;5]; % a column vector representing our set of
				    % unknowns. Used to track any row swaps.

disp('starting with the Matrix:')
A
for i = 1:size(A,[1]) % A is number of rows
	sprintf('on column %d', i)
	dgnl = A(i,i);
	n = i+1 % to be used for swapping 0 on diagonal
	while dgnl == 0 % diag is our pivot. If it's zero, elim. won't work
		A([i n],:) = A([n i],:)
		x([i n],:) = x([n i],:)
		b([i n],:) = b([n i],:)
		dgnl = A(i,i);
		n = n + 1;
	end
	for j = i+1:size(A,[2])-1
		factor = A(j,i)/dgnl;
		A(j,:) = A(j,:) - factor*A(i,:);

	end
	A
end

disp('final A in row-reduced form')
A

last_val = A(5,6)/A(5,5) % not generalized to shape !!!!
answers_reversed = [last_val]
for r = size(A,[1])-1:-1:1
	%disp('debug')
	%sprintf('working on %d', r)
	%length(answers_reversed)
	sol = A(r,end); % start with the final val of augmented matrix
	for k = 1:length(answers_reversed)
		disp('debug')
		sprintf('k is %d',k)
		sprintf('r is %d',r)
		sol = sol - A(r,size(A,[1])-k+1)*answers_reversed(k); %subtract known vals
	end
	sol = sol / A(r,r);
	answers_reversed = [answers_reversed sol]
	
end

answers = fliplr(answers_reversed)
for item = 1:length(x)
	sprintf('x%d computed as %d',x(item),answers(item))
end
disp('linsolve agrees (order has changed)')
linsolve(A0,b0)
% That was ridiculous. Would have been better in python
