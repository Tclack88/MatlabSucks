% Regression. Straight line and parabola fit
% The following data represents a stress (a) applied in N cycles until failure
N = [1 10 100 1000 10000 100000 1000000]
a = [1120 1048 991 809 643 560 425]


% Linear fit - least squares
A=[sum(N.^2) sum(N);
    sum(N) length(N)];
C=[sum(N.*a);
    sum(a)];
vec=inv(A)*C;
a1=vec(1);
b1=vec(2);
x=linspace(0,1000000);
y=a1*x+b1;
hold off
disp('plot 1')
size(N)
size(a)
plot(N,a,'ko','MarkerSize',5,'MarkerFaceColor','r')
disp('plot 2')
size(x)
size(y)
hold on
plot(x,y,'b-','LineWidth',2);
xlabel('Number of Cycles');
ylabel('Stress');

% Parabolic fit
n=length(N);
A=[n sum(N) sum(N.^2) ;
   sum(N) sum(N.^2) sum(N.^3);
   sum(N.^2) sum(N.^3) sum(N.^4)];
C=[sum(a);
    sum(N.*a);
    sum(N.^2.*a)];

figure(2)
a1=linsolve(A,C);
x=linspace(0,1000000);
y=a1(1)+a1(2)*x+a1(3)*x.^2;
hold off
plot(N,a,'ko','MarkerSize',5,'MarkerFaceColor','r')
hold on
plot(x,y,'b-','LineWidth',2);
xlabel('Number of Cycles');
ylabel('Stress');


% Neither of these are a good fit because most of the points are concentrated
% in the far left side. The scale should certainly be logarithmic
