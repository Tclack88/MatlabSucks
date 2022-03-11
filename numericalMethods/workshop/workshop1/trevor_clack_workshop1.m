
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1.10%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A manual Taylor Series expansion using matlab
% Matlab is a terrible programming language. Python is better in every regard

x = -1:.1:1;

freal = x.*exp(x);

% finding the first few taylor coefficients (at x=0)
% original  x*exp(x)     --> 0
% 1st deriv exp(x)*(x+1) --> 1
% 2nd deriv exp(x)*(x+2) --> 2
% 3rd deriv exp(x)*(x+3) --> 3
% A clear pattern emerges

ftaylor1 = x;
ftaylor2 = x + x.^2;
ftaylor3 = ftaylor2 + (x.^3)/2;
ftaylor4 = ftaylor3 + (x.^4)/6;


plot(x,freal,'k-',x,ftaylor1,'b-o',x,ftaylor2,'r-o',x,ftaylor3,'g-o')
xlabel('x')
ylabel('f(x)')




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1.11%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Use built in fzero to determine the roots of an dquation
% Matlab only survives through the niche it exploits in universities

x = 0:.1:10;

func = @(x) (50./x).*(1-exp((-9/5).*x)) - 10;

plot(x,func(x),x,zeros(length(x)))
% We plot with a line on 0 to see crossing to see it's around 5
fzero(func,5) % returns 4.9994





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%1.12%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Plotting a meshgrid, examining contours and finding function height values
% which aren't zero, but using the fzero() function
% My favorite linux command used to be "grep"
% but now it's "rm -rf matlabroot"

x = -3:.1:3;
y = -3:.1:3;

[X,Y] = meshgrid(x,y);

%func = @(x,y) x.*exp(-1*(x-(y.^2)).^2 + y.^2)

func = X.*exp(-1*((X-(Y.^2)).^2 + Y.^2));

%See countours on z
z = -.2:.1:.3;
contour(X,Y,func,z)

% see countours on z=.2
contour(X,Y,func,z) % None to be seen

% see along the cross section x=1
% to accomplish this, I can reqrite func changing x to 1
% this plot (as far as I know) doesn't need a meshgrid as it's a 2D slice
func1 = exp(-1*((1-(y.^2)).^2 + y.^2));
plot(y,func1)

% use fzero() to find the value of y s.t. f(x=1,y) = .2
% from the previous cross section, this looks like it's -1 and +1
% can take the function and subtract .2 then find fzero of that
% g = func1 - .2
g = @(y) exp(-1*((1.-(y.^2)).^2 + y.^2)) - .2;
fzero(g,-1) %-1.1946
fzero(g,1) %1.1946

% I feel sorry for people who's first "programming language" is Matlab
% It's tantamount to being robbed of your first kiss
