% Regression. Straight line fit but on a logarithmic scale
% The following data represents a stress (a) applied in N cycles until failure
N = [1 10 100 1000 10000 100000 1000000]
a = [1120 1048 991 809 643 560 425]
N = log10(N)
a = log10(a)

% Linear fit on log scale- least squares
A=[sum(N.^2) sum(N);
    sum(N) length(N)];
C=[sum(N.*a);
    sum(a)];
vec=inv(A)*C;
a1=vec(1);
b1=vec(2);

x=logspace(0,6);
y=a1*x+b1;
hold off
loglog(N,a,'ko','MarkerSize',10,'MarkerFaceColor','r')
hold on
loglog(x,y,'b-','LineWidth',2);
xlabel('Number of Cycles');
ylabel('Stress');

% This is a much better fit
vec
