A = [5 2 1; 10 -4 3; -4 9 -50];
c = [5;9;-13];

% gauss-seidel
disp('gauss-seidel')
[x1 n1] = gauss_seidel(A,c);
fprintf('After %d iterations, found value of x as:', n1)
x1

% point jacobi
disp('point jacobi')
[x2 n2] = point_jacobi(A,c);
fprintf('After %d iterations, found value of x as:', n2)
x2
disp('i.e. did not converge')

% LU decomposition
[L,U]=LUDecomposition(A);

Lc = [L c]; %L Augmented with c to do back substitution

back_sol = back_sub(L,c);
forward_sol = forward_sub(U, back_sol)

disp('ALL 3 methods converged to a solution')

%%%% functions %%%%
function [x n] = gauss_seidel(A,c)
        M = tril(A);
        N = M - A;
        x = [0;0;0]; %initial guess
        change = 100;
        n = 0 ;
        while change > 1e-6
                x_old = x;
                x = inv(M)*N*x + inv(M)*c;
                change = abs(x - x_old);
                n = n + 1;
        end

end



function [x n] = point_jacobi(A,c)
        M = diag(diag(A));
        N = M - A;
        x = [0;0;0]; %initial guess
        change = 100;
        n = 0 ;
        while change > 1e-6
                x_old = x;
                x = inv(M)*N*x + inv(M)*c;
                change = abs(x - x_old);
                n = n + 1;
        end

end


function [L,U]=LUDecomposition(A)
	% I wrote this for a previous workshop/hwk problem
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


function answer = back_sub(L,c)
	% I wrote this for a previous workshop/hwk problem
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
	% I wrote this for a previous workshop/hwk problem
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


