% Applying Root finding to the cone of a sonic boom
% the equation relating wedge angle o to oblique angle B, 
% mach number M and specific heat ratio a is given by:
% tan(o) = 2cot(B)* (M^2 (sinB)^2) - 1) / (M^2(a + cos(2B)) + 2)

% Let's plot f(B) vs B for the following values:
M = 1.4;
o = [4,8,12]; % degrees)

B = linspace(asind(1/M),90); % a range of B for plotting

hold on
for i=1:length(o)
	f_out = f(B,M,o(i));
	plot(B, f_out)
end
hold off
title('roots (zeros) of sonic boom cone by wedge angle at Mach 1.4')
legend('4 degrees', '8 degrees', '12 degrees')



% Now let's plot f(B) vs B for these next values:
M = 4;
o = [20,35,50]; % degrees)
figure
hold on
for i=1:length(o)
	f_out = f(B,M,o(i));
	plot(B, f_out);
end
hold off
title('roots (zeros) of sonic boom cone by wedge angle at Mach 4')
legend('20 degrees', '35 degrees', '50 degrees')


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
