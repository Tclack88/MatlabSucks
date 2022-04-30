% Verify Theta - Beta - M plots with contour plots

M_range = [1.4 3 4 5];
o = linspace(0,45); % A slew of ranges to examine and plot all overlaps

tiledlayout(2,2)
for M = M_range
	B_range = linspace(asind(1/M), 90);

	[X,Y] = meshgrid(B_range,o);
	nexttile
	contour(X, Y, f(X,M,Y),[0,0]);
	title('Beta plotted against Theta for Mach number',M)
	xlabel('beta')
	ylabel('theta')
end

disp('graphically we can estimate where the max theta occurs for each mach number:')
disp('M=1.4, theta max ~ 9.5')
disp('M=3, theta max ~ 34.1')
disp('M=4, theta max ~ 38.6')
disp('M=5, theta max ~ 40.9')
disp('these closely match the results obtained previously')


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
