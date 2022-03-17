% sea level will rise to a height of 35 after it rains
% H(x,y) = (y+6)^2 + 4x^2 -x^2y  by way of a contour plot
x = -3:.01:3;
y = -3:.01:3;

[X,Y] = meshgrid(x,y);

H = (Y+6).^2 + 4*X.^2 - Y.*X.^2;
z = 10:10:100;
zz = [-1000 35 1000] % other bounds given out of scale to make array
contour(X,Y,H,z) 
hold on
contour(X,Y,H,zz, 'ShowText', 'on', 'LineWidth', 3)
hold off
xlabel('x');
ylabel('y');
xlim([-3 3]);
ylim([-3 3]);
title("sea level as a function of position. 35 emphasized")
