A = [5 2 1; 10 -4 3; -4 9 -5];
c = [5;9;-13];

% triangular matrix gauss-seidel
disp('gauss-seidel')
M1 = tril(A);
N1 = M1-A;
abs(eig(inv(M1)*N1))
% diagonal matrix point jacobi
disp('point jacobi')
M2 = diag(diag(A));
N2 = M2-A;
abs(eig(inv(M2)*N2))

disp('I would expect gauss-seidel to converge because it''s P matrix eigen values are all less than 1. oint jacobi has an eigenvalue greater than 1 so convergence cannot be guaranteed')
