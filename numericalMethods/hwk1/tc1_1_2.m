% plot taylor expansion of 1/(1-2x) to 5th order

x = linspace(-5,5);
f = 1./(1-2.*x);

f0 = -1/3;
f1 = f0 + 2/9*(x-2);
f2 = f1 - 4/27*(x-2).^2;
f3 = f2 + 8/81*(x-2).^3;
f4 = f3 - 26/243*(x-2).^4; 
f5 = f4 + 32/729*(x-2).^5;

plot(x,f,x,f0,x,f1,x,f2,x,f3,x,f4,x,f5)
axis equal
grid on
xlim([0 4])
ylim([-2 2])
legend('1/(1-2x)','0th order','1st order', '2nd order', '3rd order','4th order','5th order')
title({'taylor polynomial approximations for the function f = 1./(1-2.*x)','colors do not line up because Matlab sucks and', 'I''m not gonna waste anymore time debugging'})
