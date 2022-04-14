% In this code, we will approximate the error function using gauss-legendre
% quadrature and comparing it to Matlab's in-built erf

% A provided piece of code gives the gauss-legendre points.
% It provides k points between -a and a and outputs the weights

a = 1
k = 5

ef = error_function(a,k); % "my error function" (though a big part was given)
actual = erf(a); % matlab built in
% Let's calculate the percent difference from our approximation and the builtin
percent_diff = 100*abs(ef-actual)/actual;
sprintf("the percent difference between our approximation and the built-in is %f %",percent_diff)


% Error function is used in applied mathematics such as determining the
% flow created by an infinitely long plate accelerating instantly from 
% rest to a velocity of U given by:
% u = U(1-erf(y/sqrt(4vt)))   (where v =viscosity, U=final velocity)
% Let's plot within domain (y,t) = [0,2],[0,5] with v=1 and U =1

%y = linspace(0,2);
%t = linspace(0,5);
%[Y,T] = meshgrid(y,t);
%error_function(Y./sqrt(4.*T)) 

% y = linspace(0,2,100)
% t = linspace(0,5,100)
% 
% U=1;
% v=1;
% [Y,T] = meshgrid(y,t);
% 
% erf_input = Y./sqrt(4.*v.*T);

%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%

function ret = error_function(a,k)
	% k = number of gauss_legendre points to be computed
	% symmetric endpoints of integration (-a to +a)
	% returns our approximation to the error function
	[w,x] = gauss_legendre_points(k) % provided. returns weights (w)
					 % and quadrature points (x)
	f = @(x) exp(-x.^2) % the function whose integral is being approximated
	ret = 1/sqrt(pi)*(a).*sum(w.*f(x)) % bounds are symmetric a--a/2 = a 	
end

function ret = error_function_3d(a,k)
	% k = number of gauss_legendre points to be computed
	% symmetric endpoints of integration (-a to +a)
	% returns our approximation to the error function
	[w,x] = gauss_legendre_points(k) % provided. returns weights (w)
					 % and quadrature points (x)
	y = linspace(0,2);
	t = linspace(0,5);
	[Y,T] = meshgrid(y,t);
	(Y./sqrt(4.*T)) 
	ret = 1/sqrt(pi)*(a).*sum(w.*f(x)) % bounds are symmetric a--a/2 = a 	
end
