% Multiple Integration using Simpson's rule and a meshgrid

ax = -10; % left endpoint for x
bx = 10; % right endpoint for x
n = 101; % number of intervals
x = linspace(ax,bx,n);

ay = -15; % left endpoint for y
by = 15; % right endpoint for y
y = linspace(ay,by,n);

[X, Y] = meshgrid(x,y);

F = f(X,Y); % get height of function at every point
g = simpsonsx(F, ax, bx, n);
a = simpsonsy(g, ay, by, n);
disp('area approximation with 101 divisions')
a


ns = 5:10:315;
my_vals = zeros(1,length(ns));
for i = 1:length(ns);
	n = ns(i);
	x = linspace(ax,bx,n);
	y = linspace(ay,by,n);
	[X, Y] = meshgrid(x,y);
	F = f(X,Y); % get height of function at every point
	g = simpsonsx(F, ax, bx, n);
	a = simpsonsy(g, ay, by, n);
	my_vals(i) = a;
end
h_line = 251.2*(ones(size(ns)));
plot(ns,my_vals,ns,h_line)
legend('error from my integral','actual value','location','southeast')
title('calculated value as n increases')
xlabel('n (number of divisions)')
ylabel('calculated value')


%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%

function ret = f(X,Y)
	% the function to be evaluated by simpson's rule,
	ret = 4*exp(-((X/4).^2+(Y/5).^2));
end


function area = simpsonsx (f, a, b, n)
        % approximates area under curve using simplson's rule
        % for some function f on interval a,b with n subdivisions
        range = linspace(a,b,n);

        %delta = range(2) - range(1) % this didn't work for some reason
        delta = (b-a)/n;
        mp = range(2:end-1); %mid points
        mp_odd = 3:2:n-1;
        mp_even = 2:2:n-1;
        area = delta/3*(f(:,1) + 2*sum(f(:,mp_odd),2) + 4*sum(f(:,mp_even),2) + f(:,end));
end


function area = simpsonsy (f, a, b, n)
        % approximates area under curve using simplson's rule
        % for some function f on interval a,b with n subdivisions
        range = linspace(a,b,n);
        %delta = range(2) - range(1); % this didn't work for some reason
        delta = (b-a)/n;
        mp = range(2:end-1); %mid points
        mp_odd = 3:2:n-1;
        mp_even = 2:2:n-1;

        area = delta/3*(f(1) + 2*sum(f(mp_odd)) + 4*sum(f(mp_even)) + f(end));
end
