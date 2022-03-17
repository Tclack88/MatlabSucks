% Height in a cylindrical water tank laid on its side.
% Volume of water 100m^3m, L=50m
% Matlab blows

% Volume is given by: (r^2*acos((r-h)/r) - (r-h)*sqrt(2*4*h-h^2))*L
% for a radius r, 0 < h < 2r. rmin = sqrt(V/(Lpi))
% substituting L = 50 and knowing we are interested in V = 100,
% the roots we are trying to find are of
% f(h) = (r^2*acos((r-h)/r) - (r-h)*sqrt(2*4*h-h^2))*50 - 100

% plot from 0.5 < h < 1.5 for r = 1,2,5


x = linspace(0.5,1.5);
line([.5,1.5], [0,0])
%plot(x,one_bar)
hold on
for r = [1 2 5]
	plot(x,f(r, x));
end
xlabel('water height (m)')
ylabel('volume deviation from 0 (which is 100m^3)')
title({'volume of a cylindrical pool as a function of height for different radii';'Horizontal line shows available water volume'})
legend('100m^3 mark', 'r=1', 'r=2','r=5')
hold off

% approximate this height for r=2 using method of False Position
% root at r = 2 graphically appears somewhere between .7 and 1

disp('use method of false position to find root at r=2, i.e. numerically calculate at r=2 the height of a cylindrical swimming pool')
false_position(.7,1,@f,2, x);

% Now we very r and see what the ranges of heights are
clear
figure
x = linspace(1.4,10); % new set of values
line([1.4,10], [0,0])
%plot(x,one_bar)
hold on
for h = [.5 .7 1]
	plot(x,f(x, h));
end
xlabel('radius of tank (m)')
ylabel('volume deviation from 0 (which is 100m^3)')
title({'volume of a cylindrical pool as a function of tank radius for different water height';'Horizontal line shows available water volume'})
legend('100m^3 mark', 'h=.5', 'r=.7','r=1')
hold off

% at r=.7, h is somewhere between 3 and 4.
% Let's check with the method of false position
disp('use method of false position to find root at h=.7, i.e. numerically calculate at h=.7 the radius of a cylindrical swimming pool')
false_position(3, 4, @f, x, .7) 

function result = f(r, h)
	result =  (r.^2.*acos((r-h)./r) - (r-h).*sqrt(2.*r.*h-h.^2))*50 - 100;
end

function res = false_position(x_l, x_u, f, r, h)
	% var is r or h, depending on the function inputted
	eps = 1.0e-6;
	n = 0;
	if size(r) == 1
		var = r;
		sprintf('holding radius constant at %d',var)
		f_l = f(var, x_l);
		f_u = f(var, x_u);
		x_r = x_u - f_u*(x_l-x_u)/(f_l - f_u); % get initial x_r
		while abs(f(var, x_r))>eps
		    f_l=f(var, x_l);
		    f_u=f(var, x_u);
		    f_r=f(var, x_r);
		    if f_r*f_l < 0
		        x_u= x_r;
		    elseif f_r*f_u < 0
		        x_l= x_r;
		    end
		    x_r = x_u - f_u*(x_l-x_u)/(f_l - f_u);
		    n = n + 1;
		end
	elseif size(h) ==1
		var = h;
		sprintf('holding height constant at %d',var)
		f_l = f(x_l, var);
		f_u = f(x_u, var);
		x_r = x_u - f_u*(x_l-x_u)/(f_l - f_u); % get initial x_r
		while abs(f(x_r, var))>eps
		    f_l=f(x_l, var);
		    f_u=f(x_u, var);
		    f_r=f(x_r, var);
		    if f_r*f_l < 0
		        x_u= x_r;
		    elseif f_r*f_u < 0
		        x_l= x_r;
		    end
		    x_r = x_u - f_u*(x_l-x_u)/(f_l - f_u)
		    n = n + 1;
		end
	else
		disp('error, one of r or h must be a single number')
	end
sprintf("method of false position finds root of %f after %d iterations", x_r,n)
res = x_r;
end
