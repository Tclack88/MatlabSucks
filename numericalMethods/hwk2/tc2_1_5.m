% We found that at M = 4, omax occurs around o = 66. Let's now try to find
% Bl and Bu as we keep M =4 but vary theta o from 0 to omax (approx 37)

M = 4
o = linspace(0,39); % saw 37 as initial approximation where o was zero. Here
% I have slowly dialed it up to find a closer value where they meet

% we found for Mach 4, when o = 0, BL = arcsin(1/4), Bu = 90. Then we saw
% graphically that as o grew, the f(B) plot moved down. Bringing BL and Bu 
% closer. As such, we can use this initial value as our lower and upper range

B_range = linspace(asind(1/M), 90);

% We also can see that the peak of the plot occurs around B = 68, we can use
% a bisection method to find the roots then given 0 and 68 as the initial
% l and u bounds for the left root, and 68 and 90 for the initial l and u bounds
% for the right root.
Ll = asind(1/M);
Lu = 68;
Rl = 68;
Ru = 90;

BLs = []; % an array to store all the BLs
BUs = []; % an array to store all the BUs

for i = 2:length(o)-1;
	BL = bisection(@f,M,o(i), Ll, Lu);
	BU = bisection(@f,M,o(i), Rl, Ru);
	BLs(end+1) = BL;
	BUs(end+1) = BU;
end

plot(o(2:end-1),BLs, o(2:end-1),BUs)
xlabel('theta (degrees)')
ylabel('Beta (degrees)')
title('The value for BL and BU for different thetas')
legend('BL', 'BU')




%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%

function ret = f(B,M,o)
	% B - oblique angle
	% M - Mach number
	% o - wedge angle
	a = 1.4; % ratio of specific heats
	num = M^2 .* sind(B).^2 - 1;
	denom = M^2 * (a + cosd(2.*B)) + 2;
	ret = 2 .* cotd(B).*num./denom - tand(o);
end


function x_r = bisection(f, M, o, x_l, x_u)
	eps = 1.0e-3; % Accuracy threshold
	maxiter = 50; % maximum number of iterations
	
	x_r=(x_l+x_u)/2;
	n = 0;
	f(x_r, M, o)
	while abs(f(x_r, M, o))>eps
	    f_l=f(x_l, M, o);
	    f_u=f(x_u, M, o);
	    f_r=f(x_r, M, o);
	    if f_r*f_l < 0
	        x_u= x_r;
	    elseif f_r*f_u < 0
	        x_l= x_r;
	    end
	    x_r=(x_l+x_u)/2;
	    x_r
	    f(x_r,M,o)
	    n = n + 1; % augment to count
	    if n == maxiter
		    break
	    end
	end
end

