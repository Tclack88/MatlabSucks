% Multiple Integration using Simpson's rule and a meshgrid

ax = 0; % left endpoint for x
bx = 5; % right endpoint for x
n = 101; % number of intervals
x = linspace(ax,bx,n);

ay = -3; % left endpoint for y
by = 3; % right endpoint for y
y = linspace(ay,by,n);

[X, Y] = meshgrid(x,y);

F = f(X,Y); % get height of function at every point
g = simpsonsx(F, ax, bx, n);
a = simpsonsy(g, ay, by, n);
disp('area approximation with 101 divison points')
a

% compare as the number of divisions increases
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
h_line = 1456*(ones(size(ns)));
plot(ns,my_vals,ns,h_line)
legend('error from my integral','actual value', 'location','southeast')
title('calculated value as n increases')
xlabel('n (number of divisions)')
ylabel('calculated value')


%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%

function ret = f(X,Y)
	% the function to be evaluated by simpson's rule,
	ret = X.^2 - 2*Y.^3 + X.*Y.^4;
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

function int = Simpson(f,x)
	delta = x(2)-x(1);
	delta = (x(end) - x(1)) / length(x)
	int = f(1)
	int = int + int+4*sum(f(2:2:end-1));
	int = int + int+2*sum(f(3:2:end-1));
	int = int + f(end);
	int = delta/3*int;
end
