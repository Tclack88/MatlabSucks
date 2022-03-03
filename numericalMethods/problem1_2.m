% Use built in fzero to determine the roots of an dquation
% Matlab only survives through the niche it exploits in universities

x = 0:.1:10;

func = @(x) (50./x).*(1-exp((-9/5).*x)) - 10;

plot(x,func(x),x,zeros(length(x)))
% We plot with a line on 0 to see crossing to see it's around 5
fzero(func,5) % returns 4.9994
