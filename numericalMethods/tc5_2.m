clear
N = 40;
x = linspace(0,20,N);
dx = x(2) - x(1)

dts = [.01 .05 .1 .2 .5 1];
figure
hold on
for dt = dts;
        O = initiate_phis(x);
        evals = get_eig_vals(O,dx)';
        X = real(evals)*dt;
        Y = imag(evals)*dt;
        plot(X,Y,'o');
end
legend_strings = "dt: " + string(dts);
legend(legend_strings,'location','southwest')
title('stablility plots for various dt time steps with RK2')



relamdt=-3:0.1:1;
imlamdt=-2:0.1:2;
[RElamdt,IMlamdt]=meshgrid(relamdt,imlamdt);
axis square;
lamdt=RElamdt+i*IMlamdt;
sig=1+lamdt+lamdt.^2/2;
contour_levels=[1 0 0 0];
contour(RElamdt,IMlamdt,abs(sig),contour_levels,'linewidth',4);
xlabel('Re \lambda \Delta t');
ylabel('Im \lambda \Delta t');
hold on;

