% at M = 4, when o=20, f(B) = 0 around B = 30 and B = 85.
% But we can check with a root finding method. I'll do the bisection method 
% because I don't need to find derivatives or rearrange to solve for some B
% given this, I will conserviatively choose my bounds as 25-35 and 80-90

xl = 25;
xu = 35;
bisection(@f, xl, xu)

xl = 80;
xu = 90;
bisection(@f, xl, xu)

%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%

function ret = f(B)
	% B - oblique angle
	% o - wedge angle
	M = 4; % Mach number
	o = 20; % wedge angle - degrees
	a = 1.4; % ratio of specific heats
	num = M^2 * sind(B)^2 - 1;
	denom = M^2 * (a + cosd(2*B)) + 2;
	ret = 2 * cotd(B)*num/denom - tand(o);
end

function  ret = bisection(f, xl, xu)
	gap = 10;
	n = 0;
	while gap > .0000001
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

