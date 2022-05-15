% Initial value problem. Sailboat in wavy waters with Runga Kutta
% function is 100 x'' = 50Vcos(o) + 6Asin(x/10) - 60x'




for n=1:N
	k1=f(n);
	k2=(f(n)+(h/2)*k1);
	f(n+1)=f(n)+((1/2)*k1+(1/2)*k2)*h;
	t(n+1)=t(n)+h;
end;


% The equation can be rearranged as:
% x'' = 1/2Vcos(o) + 3/50Asin(x/100) - 3/5x'
% Here, we can reduce to a system of equations:
% u = x  u'=v
% v = x' v'= 1/2Vcos(o) + 3/50Asin(u/100) - 3/5v
%   0          1
% -3/5   1/2Vcos(o) + 3/50Asin(u/100) 


