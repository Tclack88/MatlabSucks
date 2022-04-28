% Theta - Beta - M plot
% Previously we determined roughly what theta max is for M = 4. Let's
% repeat for other values of M (1.4, 2, and 5).

M_range = [1.4 3 4 5];
o = linspace(0,45); % A slew of ranges to examine and plot all overlaps


tiledlayout(2,2)
for M = M_range
	B_range = linspace(asind(1/M), 90);

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

	nexttile
	plot(o(2:length(BLs)+1),BLs, o(2:length(BUs)+1),BUs)
	xlabel('theta (degrees)')
	ylabel('Beta (degrees)')
	title('The value for BL and BU for different thetas for M=',M)
	legend('BL', 'BU')
end


disp('graphically we can see theta max occurs at the following values of M:')
disp('M=1.4 -> theta = 9')
disp('M=3 -> theta =33.6')
disp('M=4 -> theta =38.6')
disp('M=5 -> theta =40.9')



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
	eps = 1.0e-6; % Accuracy threshold
	maxiter = 50; % maximum number of iterations
	
	x_r=(x_l+x_u)/2;
	n = 0;
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
	    n = n + 1; % augment to count
	    if n == maxiter
		    break
	    end
	end
end

