% using the built-in fzero to verify the manual results obtained previously
gx = @(x) sin(2*x) + x^2/3 - 3;

fzero(gx, 3.1)
