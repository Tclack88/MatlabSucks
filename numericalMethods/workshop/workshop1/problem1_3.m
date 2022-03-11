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
