% Multiple Integration using Simpson's rule and a meshgrid

ax = 0; % left endpoint for x
bx = 5; % right endpoint for x
n = 8; % number of intervals
x = linspace(ax,bx,n);

ay = -3; % left endpoint for y
by = 3; % right endpoint for y
y = linspace(ay,by,n);

[X, Y] = meshgrid(x,y);

F = f(X,Y); % get height of function at every point
size(F)

%inner_integral = simpsons(F, ax, bx, n) % integrated over x

g = (bx-ax)/n * trapz(F)
%(by-ay)/n *trapz(g)
%trapz(x,F,2)
%trapz(y, trapz(x,F,2))



%my_area = simpsons(@f, a, b, n)


%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%

function ret = f(X,Y)
	% the function to be evaluated by simpson's rule,
	ret = X.^2 - 2*Y.^2 + X.*Y.^4;
end


% function area = simpsons (f, a, b, n)
%         % approximates area under curve using simplson's rule
% 	% for some function f on interval a,b with n subdivisions
%         range = linspace(a,b,n);
%         delta = range(2) - range(1);
%         mp = range(2:end-1); %mid points
%         mp_odd = 1:2:length(mp);
% 
%         mp_even = 2:2:length(mp);
% %	mp_odd
% %	mp_even
% f(1)
% 
% 	f(mp_odd(1,:))
% 
% 
%         area = delta/3*(f(1) + 4*sum(f(mp_odd,:)) + 2*sum(f(mp_even,:)) + f(end));
% end
