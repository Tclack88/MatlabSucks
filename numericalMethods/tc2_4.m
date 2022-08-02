k = 100

for n = [3, 100, 10000];
	n
	iterative(n,k)
	LU(n)
end




function res = iterative(n,k)
	res = (n^2)*k;
end

function res = LU(n)
	res = (n^3);
end
