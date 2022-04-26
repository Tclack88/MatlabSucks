% Integration using Simpson's rule

a = 1; % left endpoint
b = 10; % right endpoint
n = 100; % interval size
X = linspace(a,b,n); % 10 subdivisions (arbitrary choice)


area = .7122645; % A provided "actual" area
matlab_area = trapz(X, f(X)) % the "actual value" from a matlab built-in
my_area = simpsons(@f, a, b, n)

disp('percent error from my function:')
100*abs(my_area - area)/area

disp('percent error from matlab function:')
100*abs(matlab_area - area)/area

disp('clearly matlab''s built-in is more accurate (with this chosen interval')


%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%

function ret = f(x)
	% the function to be evaluated by simpson's rule
	ret = sin(x)./x;
end


function area = simpsons (f, a, b, n)
        % approximates area under curve using simplson's rule
	% for some function f on interval a,b with n subdivisions
        range = linspace(a,b,n);
        delta = range(2) - range(1);
        mp = range(2:end-1); %mid points
        mp_odd = mp(1:2:end);

        mp_even = mp(2:2:end);

        area = delta/3*(f(a) + 4*sum(f(mp_odd)) + 2*sum(f(mp_even)) + f(b));
end



