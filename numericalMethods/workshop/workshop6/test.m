xi=[1 4 6];
fxi=[2 3.386294 3.791760];
x=2
N=length(xi);
fx=0;
for i=1:N
Lix=1.0;
for j=1:N
if j ~= i
Lix=Lix*(x-xi(j))/(xi(i)-xi(j));
end
end
fx=fx+fxi(i)*Lix;
end
