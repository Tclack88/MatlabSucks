%LU decomposition
% system of linear equations:
% 4x - y - 2z = 1
% x - 2y + 3z =4
% x + y +z = 1

A = [4 -1 -2;
	1 -2 3;
	1 1 1];

c = [1; 4; 1];

% Check solutions:
out = A\c

% Decompose A into LU using crout's Method


[L,U]=LUDecomposition(A)
% Here is A broken into an L and U
[A L U]

Lc = [L c]; %L Augmented with c to do back substitution

R= LowerTriangularSolver(L,c);
x = UpperTriangularSolver(U,R)

% % Back sub adapted from forward sub I did with gaussian elimination
% % in workshop3
% back_sol = back_sub(L,c)
% % This agrees with the built-in tool
% L\c
% 
% Uc = [U back_sol]; % augmented matrix to do back substitution
% 
% % Here's the forward sub I tried for a guassian elimination in workshop 3
% forward_sol = forward_sub(U,back_sol)
% % This agrees with the built-in tools
confirmation = U\(L\c)
% 
% 
% % test with a different solution set
c2 = [2;6;7];
% 
% back_sol2 = back_sub(L, c2);
% 
% forward_sol2 = forward_sub(U, back_sol2)
R2= LowerTriangularSolver(L,c2);
x2 = UpperTriangularSolver(U,R2)
confirmation = U\(L\c2) %check output with built-in tool

%%%%%%%%%%%%%%%% back and forward substitution functions %%%%%%%%%%%%%%
function answer = back_sub(L,c)
	Lc = [L c];
	m = size(L,[1]);

	first_val = Lc(1,m+1)/Lc(1,1);
	answers = [first_val];
	for r = 2:m
        	val = Lc(r,end); % start with the final val of augm. matrix
        	for k = 1:length(answers)
                	val = val - Lc(r,k)*answers(k); %subtract known vals
        	end
        	val = val / Lc(r,r);
        	answers = [answers val];

	end
	answer = answers';
end


function answer=forward_sub(U,c)
	Uc = [U c]; % augmented matrix
	m = size(U,[1]); % side length of square matrix

	last_val = Uc(m,m+1)/Uc(m,m);
	answers_reversed = [last_val];
	for r = m-1:-1:1
        	val = Uc(r,end); % start with the final val of augm. matrix
        	for k = 1:length(answers_reversed)
                val = val - Uc(r,m-k+1)*answers_reversed(k); %subtract known vals
        	end
        	val = val / Uc(r,r);
        	answers_reversed = [answers_reversed val];
	end
	answer = fliplr(answers_reversed)';
end

function [L,U]=LUDecomposition(A)
L = zeros(size(A));
U = eye(size(A));

[m n] = size(A);

L(:,1) = A(:,1); % first col  of L is first col of A

for j=2:n % First row of U is 1st row of A elements/ L(1,1)
	U(1,j)=A(1,j)/L(1,1);
end

for j=2:n
	for i=j:n % get the L's which is A-(sum LU products)/(U in diag=(1))
		sum = A(i,j);
		for k=1:j-1
			sum = sum - L(i,k)*U(k,j);
		end;
		L(i,j) = sum;
	end;
	for k=j+1:n % get U's which is A - (sum LU products)/L
		sum = A(j,k);
		for i=1:j-1
			sum = sum - L(j,i)*U(i,k);
		end;
		U(j,k)=sum/L(j,j);
	end;
end;
end

function x=UpperTriangularSolver(A,C)
[n,n]=size(A); % Finding the size of the problem
x=zeros(n,1); % set up  x vector with elements initially set to zero
x(n)=C(n)/A(n,n);
 for i=n-1:-1:1
    sum=0;
    for j=i+1:n
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(C(i)-sum)/A(i,i);
end
end


function x=LowerTriangularSolver(A,C)
[n,n]=size(A); % Finding the size of the problem
x=zeros(n,1); % set up  x vector with elements initially set to zero
x(1)=C(1)/A(1,1);
for i=2:n
    sum=0;
    for j=1:i-1
        sum=sum+A(i,j)*x(j);
    end
    x(i)=(C(i)-sum)/A(i,i);
end
end
