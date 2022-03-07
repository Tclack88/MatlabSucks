
%f(x) = sin(2x) + x^2/3 - 4
%plot from -5 to 5

x = linspace(-5,5);
fx = sin(2*x) + x.^2/3 - 4;
plot(x,fx)

hold on
scatter(-2.5,-1,10,'o','r','MarkerFaceColor', 'r') % guess for lower f(x)=-1
scatter(3.1,-1,10,'o','r','MarkerFaceColor', 'r')  % guess for upper f(x)=-1
hold off
grid on
set(gca,'xtick',[-5:1:5]);

