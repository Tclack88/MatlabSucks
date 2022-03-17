% method of bisection compared to false position

% bisection method:
disp('bisection')
f=@(x) x.^3-x;

x_l = -2.2; %Initial lower guess of root
x_u = -0.8; %Initial upper guess of root


eps = 1.0e-6; %Accuracy of your answer
maxiter = 50; %maximum number of iterations

x_r=(x_l+x_u)/2

while abs(f(x_r))>eps
    f_l=f(x_l);
    f_u=f(x_u);
    f_r=f(x_r);
    if f_r*f_l < 0
        x_u= x_r;
    elseif f_r*f_u < 0
        x_l= x_r;
    end
    x_r=(x_l+x_u)/2
end

% false position
disp('false position')

x_l = -2.2; %Initial lower guess of root
x_u = -0.8; %Initial upper guess of root


eps = 1.0e-6; %Accuracy of your answer
maxiter = 50; %maximum number of iterations

x_r = 80;

while abs(f(x_r))>eps
    f_l=f(x_l);
    f_u=f(x_u);
    f_r=f(x_r);
    if f_r*f_l < 0
        x_u= x_r;
    elseif f_r*f_u < 0
        x_l= x_r;
    end
    x_r=fp_xr(x_l, x_u, f_l, x_u)
end


% matlab requires functions to be at the bottom. This is not a limitation
% in python, which sensibly just requires a function to be defined before
% being called. From top to bottom, linearly, intuitively.
function x_r = fp_xr(x_l, x_u, f_l, f_u)
       x_r = x_u - f_u*(x_l-x_u)/(f_l - f_u)
end
