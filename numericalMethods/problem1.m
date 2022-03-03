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



