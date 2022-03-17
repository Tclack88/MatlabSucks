% Water shear and reynold's number
% The colebrook equation for a smooth pipe is:
% g(f, Re) = 1/sqrt(f) + 2*log10(2.51/(Re*sqrt(f))) = 0
% we wish to numerically find f for different values of Re

f_range = linspace(0, .1);
line([0 .1], [0 0]);
hold on
for Re = [10e4, 10e5, 10e6]
	plot(f_range, g(f_range,Re,0))
end
%xlim([.01 .1])
title('colebrook output as a function of Reynolds number')
xlabel('friction factor')
ylabel('g')

% The haaland equation is used to estimate f and is given by:
% 1/sqrt(f) ~ -1.8*log10(6.9/Re] (for a smooth pipe)
% as such, f ~ (-1.8*log10(6.9/Re])^-2 
% Let's plot for the same set of Reynold's numbers:
f_factors = []; % save for later (when doing fixed point iteration)
for Re = [10e4, 10e5, 10e6]
	f_val = f(Re, 0);
	f_factors = [f_factors f_val];
	sprintf("For Re= %d, f= %d by haaland equation", Re, f_val)
	xline(f_val);
end

% These lines match up very nicely proving the haaland equation
% to be a good approximation
% Let's use these initial guesses to find the actual value
% since this is a single guess, I will use fixed-point iteration


Re_num = [10e4 10e5 10e6]
for k = 1:3
	fac = f_factors(k);
	Re = Re_num(k);
	fixed_point_iter(@g, @fp_g, Re, 0, fac, 1e-10);
end

% Let's do this with more values of Re and plot f as a function of Re
Re_range = logspace(3, 8);
f_range_pred = f(Re_range, 0);
f_range_calc = [];
for k = 1:length(Re_range);
	Re = Re_range(k);
	fac = f_range_pred(k);
	calc_f = fixed_point_iter(@g, @fp_g, Re, 0, fac, 1e-10);
	f_range_calc = [f_range_calc calc_f];
end
figure;
% plot on a log-scaled plot for the x values
loglog(Re_range, f_range_calc);
title('f factors as a function of Re')
xlabel('Re')
ylabel('f (calculated)')


% Let's now change it so that we are not dealing with non-smooth pipes
% i.e. several non-zero Ed values
eD_range = [.00002 .01 .02 .05]
figure
%line([10e3 10e8], [.7 .7]);
for eD = eD_range
	f_range_calc_witheD = [];
	for k = 1:length(Re_range);
		Re = Re_range(k);
		fac = f(Re, eD);
		calc_f = fixed_point_iter(@g, @fp_g, Re, eD, fac, 1e-10);
		f_range_calc_witheD = [f_range_calc_witheD calc_f];
	end
	loglog(Re_range, f_range_calc_witheD);
	hold on
end
hold off
title('f factors as a function of Re (with eD factor included')
xlabel('Re')
ylabel('f (calculated)')
legend('eD=.00002', 'eD=.01', 'eD=.02', 'eD=.05')

%%%%%%%%%%%%%%%%%%% my functions %%%%%%%%%%%%%%%%%%%%%%%%%
function ret =  g(f, Re, eD)
	ret = 1./sqrt(f) + 2*log10(eD/3.7 + 2.51./(Re*sqrt(f)));
end

function ret =  fp_g(f, Re, eD)
	% fixed point iteration of g
	% 1./sqrt(f) + 2*log10(2.51./(Re*sqrt(f))) = 0
	% rearrange for an f (the left side is fine)
	ret = (2*log10(eD/3.7 + 2.51./(Re*sqrt(f))))^-2;
end

function ret = f(Re, eD)
       ret =  (-1.8*log10(6.9./Re + (eD/3.7).^1.11)).^-2; 
end

function ret = fixed_point_iter(g, fp_g, Re, eD, xi,lim)
	% xi -  initial guess
	% Re - reynolds number
	% lim limit
	l = 10; % large initial number
	n = 0;
	while l > lim
		xi = fp_g(xi, Re, eD); % get next iteration
	        l = abs(g(xi, Re, eD)); % compare to actual function
	        n = n + 1;
	end
	
	sprintf('Arrived at approximation of %f after %d cycles', xi, n)
	ret = xi;
end


