% visualize height above sea level as given by 
% H(x,y) = (y+6)^2 + 4x^2 -x^2y  by way of a contour plot
x = -3:.01:3;
y = -3:.01:3;

[X,Y] = meshgrid(x,y);

H = (Y+6).^2 + 4*X.^2 - Y.*X.^2;
z = 10:10:100;
contour(X,Y,H,z, 'ShowText', 'on') % add inline labels of contours
xlabel('x');
ylabel('y');
xlim([-3 3]);
ylim([-3 3]);
title('sea level as a function of position')
