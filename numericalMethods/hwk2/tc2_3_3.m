% From the previous analysis, A good plot was found over the log plot
% of:

N = [1 10 100 1000 10000 100000 1000000];
a = [1120 1048 991 809 643 560 425];

% values obtained previously from the vector outpouts of the previous
a1 = -.0712;
b1 = 3.0935;


% x=logspace(0,6);
% y=a1*x+b1;
% hold off
% loglog(N,a,'ko','MarkerSize',10,'MarkerFaceColor','r')
% hold on
% loglog(x,y,'b-','LineWidth',2);
% xlabel('Number of Cycles');
% ylabel('Stress');


%x = 10^b*(10.^(a1.*log(N)))
%f = (10^b1)*(10.^(a1*log(N)))

plot(log10(N),log10(a),'ko','MarkerSize',10,'MarkerFaceColor','r')
hold on

x=logspace(0,6);
f = 10.^(a1*x + b1)
plot(x,f,'b-','LineWidth',2);
