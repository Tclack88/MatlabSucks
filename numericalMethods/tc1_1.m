V = linspace(.002,.004);
fV = f(V);
plot(V,fV, V,zeros(size(V)));
title('plot of Pressure as a function of volume')
xlabel('Volume')
ylabel('Pressure')

disp('By inspection P=0 around V = .00295')

function res =  f(V)
	T = 230;
	P = 50000;
	b = .001863;
	R = .5183;
	a = 12.56;
	res = R*T./(V-b) - a./(V.*(V+b)*sqrt(T))-P; % subtract P to find when f is zero
end
