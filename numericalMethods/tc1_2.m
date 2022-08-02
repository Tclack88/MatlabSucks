clc;
clear all;
V = linspace(.002,.004);
fV = f(V);

% use bisection guessing left and right of .003
xl = .0029;
xu = .0031;
bisection(@f, xl, xu)
disp('we are asked how many kg a 3.5m^3 tank will hold')
f(2.5)


function res =  f(V)
	T = 230;
	P = 50000;
	b = .001863;
	R = .5183;
	a = 12.56;
	res = R*T./(V-b) - a./(V.*(V+b)*sqrt(T))-P; % subtract P to find when f is zero
end

function  ret = bisection(f, xl, xu)
        gap = 10;
        n = 0;
        while gap > .00000000001
                xr = (xu + xl)/2;
                U = f(xu);
                L = f(xl);
                R = f(xr);
                if R*L > 0
                        xl = xr;
                elseif R*U > 0
                        xu = xr;
                end
                n = n + 1;
                gap = xu - xl;
                xr;
        end
        ret = xr;
end
