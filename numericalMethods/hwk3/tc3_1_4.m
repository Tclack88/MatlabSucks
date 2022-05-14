% Plot of exact solution to ODE

xs = linspace(0,1);
ys = f(xs);
plot(xs,ys);
title('exact ODE plot of y''''+y''-6y = 15sin(12x)');


disp('numerical calculations match very well with the exact solution')
disp('even with "cruder" divisions, it does not deviate unreasonably')
disp('as one would expect from a Euler''s method plot.')
disp('if both methods are ran on top of each other, it woud not be clear that')
disp('they are distinct until zooming in to the scale where each')
disp('division differs in the thousanths place')

function y = f(x)
	y = -1.264936788*exp(-3*x) + 0.2728859139*exp(2*x) - 5/629*cos(12*x) - 125/1258*sin(12*x);
end
