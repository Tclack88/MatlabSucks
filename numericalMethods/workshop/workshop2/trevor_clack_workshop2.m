% Temperature on a square plate

x = linspace(-1,1);
y = linspace(-1,1);
[X,Y] = meshgrid(x,y);
T = 80*exp(-(X-1).^2).*exp(-3*(Y-1).^2);

z = [70 60 50 40 30 27 20 10 5];

contour(X,Y,T,z)
hold on
[c,h] = contour(X,Y,T,[0 27 80]) 	% emphasize 27 isotherm (80 and 0 given
h.LineWidth = 2;			% to provide an array, but now shown
xlim([0 1]);
ylim([0 1]);
xlabel('x')
ylabel('y')
title('a plot of isotherms. 27 emphasized')
hold off

% partial derivatives
%dTdx = -160*(x-1)*exp(-(x-1)^2)*exp(-3(x-1)^2)
%dTdy = -480*(y-1)*exp(-(x-1)^2)*exp(-3(y-1)^2)

% A spider follows the isotherm along 27C
% i. In this case, we are interested in which values (x,y) give the temperature
% 27C, so we set T(x,y) = 27. In order to use the Newton-Raphson method,
% we must look at the hyperfunction f = T(x,y) - k, where k is the isotherm
% we are interested in (in this case 27). Knowing that the spider is at x=0,
% this simplifies our equation to one of one variable, namely:
% f(y) = 80*exp(-(0-1)^2)*exp(-3*(y-1)^2) - 27 = 0
% f(y) = 80*exp(-1)*exp(-3*(y-1)^2) - 27 = 0
% f(y) = 80*exp(-3*(y-1)^2 - 1) - 27 = 0
% This can be solved analytically relatively easily, however we will solve
% it numerically for the sake of the exercise

fy = @(y) 80*exp(-3.*(y-1).^2 - 1) - 27
% ii plot
figure
plot(y, fy(y), y, zeros(length(y))) % plot. Set line at f(y)=0 to make initial guess
xlim([0 1]);
%ylim([0 1]);
xlabel('y')
ylabel('f(x=0,y)')
title('80exp(-3(y-1)^2 -1) - 27')
%Tt = 80*exp(-(0-1)^2)*exp(-3*(1.1695-1)^2)

% Graphically, it appears that f(0) (i.e. T(x,y)=27) occurs at roughly y=0.8
% Newton-Raphson method
% we need our derivative (at x=0) which was found above
Fy = @(k) 80*exp(-3*(k-1)^2 - 1) - 27
dFy = @(k) -480*(k-1)*exp(-3*(k-1)^2 - 1)

NR = @(k) k - (Fy(k))/(dFy(k))

xi = .8;
l = 10; % large initial number
n = 0; 
while l > .00001
	xi = NR(xi);
	l = abs(Fy(xi));
	n = n + 1;
end

sprintf('Arrived at approximation of %f after %d cycles', xi, n)

% Spider wants to follow along the isotherm 27C from x=0 to x=1
% I don't know how to find this programatically, but we are essentially
% solving for x,y where T=27 which is as follows:
%
% 27 = 80*exp(-3*(k-1)^2 - 1)
% log(27/80) = -(x-1)^2 -3(y-1)^2
%
% 1 =    (x-1)^2    +      (y-1)^2
%      (log(80/27)     (1/3)(log(80/27))
%
% this is an ellipse centered at 1,1 with radii
% x = sqrt(log(80/27)) ~.6017 
% y = sqrt((1/3)(log(80/27))) ~1.0422
% This looks like our plot as when x=1, y is at a minimum of 1-.6017=.398
% when y = -1, x ~-.0422 leading us to belive that when x=0, y is a little less
% than 1. This can be calculated explicitle by setting x=0 in the equation
% above on line 74 and solving explicitly for y which gives 2 points
% y= 1.169 and 0.8305. (from 1 +- sqrt( (1 - 1/log(80/27)) * (1/3)*log(80/27) ))
% the latter value matches both of our plots of y,y(0) and x,y
