k = 100;
n=3;
n
disp('iterative:')
iterative(n,k)
disp('LU decomp')
LU(n)
disp('clearly LU decomp is faster for a 3x3 matrix')





function res = iterative(n,k)
	res = (n^2)*k;
end

function res = LU(n)
	res = (n^3);
end
