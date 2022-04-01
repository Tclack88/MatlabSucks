% Convergence for newton-raphson for linear and quadratic orders
load ENGR20005_Workshop5p8.mat
% The above loads data in some 1 by x matrix (x varies) for each of the
% following methods which we used across workshops 2 and 3 to estimate ...
% something for which the actual value was 3 for the following methods
%
% data_bisection
% data_false_position
% data_fixed_point_iteration
% data_newton_raphson
% data_secant
%
% We wish to find error from each method

E_b = get_error(data_bisection, 3);
E_fp = get_error(data_false_position, 3);
E_fpi = get_error(data_fixed_point_iteration, 3);
E_nr = get_error(data_newton_raphson, 3);
E_s = get_error(data_secant, 3);




% plot Error of the (i+1)th iteration against the ith iteration error
% This process is rather inefficent compared to a better programming
% language like python
figure
plot_error(E_b);
hold on
plot_error(E_fp);
plot_error(E_fpi);
plot_error(E_nr);
plot_error(E_s);
legend('bisection', 'false position', 'fixed point iteration', 'newton raphson', 'secant')
hold off


% Newton Raphson and Secant do not have linear approaches. We wish to see
% what these are by plotting exponentiall or logarithmic approach
% that is we will plot a (Ei_1) = a(Ei)^b as log(Ei_1) = log(a) + b*log(Ei)
% Then we can use the slope and y intercepts to determine the coefficients 
% a and b

figure
p_nr = poly_plot_error(E_nr)
hold on
p_s = poly_plot_error(E_s)
hold off
legend('newton raphson', 'secent')

% newton-raphson method vals:
p_nr(1)
p_nr(2)

% secant method vals:
p_s(1)
p_s(2)


%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%
function E = get_error(data, r);
	E = double(abs(data - r));
end

function plot_error(err);
	Ei = err(1:end-1);
	Ei_1 = err(2:end);
	plot(Ei, Ei_1)
end

function p =  poly_plot_error(err);
	Ei = err(1:end-1);
	Ei_1 = err(2:end);
	p = polyval(Ei, Ei_1,1);
	plot(log(Ei), log(Ei_1));
end
