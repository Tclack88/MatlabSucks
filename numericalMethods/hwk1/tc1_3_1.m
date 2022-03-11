% plot of force as a function of angle the force is applied
% for pulling a box along a ground with friction
u = .58;
g = 10;
m = 6;

F = @(t) u*m*g./(cos(t) +u*sin(t));

t_range = linspace(0,pi/2);

plot(t_range,F(t_range))
% make it pretty
title('Force as a function of angle from 0 to π/2')
xticks([0 pi/4 pi/2]);
xticklabels({'0','π/4','π/2'})
xlabel('angle (radians)')
ylabel('Force (Newtons)')

% Minimum force is not horizontal, it's actually at some angle
% this is because some vertical component will remove some of the friction
% this is because friction is a function of Normal force 
