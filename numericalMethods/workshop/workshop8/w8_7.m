% Approximate derivates and compare using
% differencing schemes and spectral difference


f_list = {@f1,@f2,@f3,@f4};
df_list = {@df1,@df2,@df3,@df4};

% approximate derivatives using differencing schemes for 10 points
n = 10;
x = linspace(-1,1,n);
err = [];
for fi = 1:length(f_list);
	y = f_list{fi}(x);
	dy = [FDS(x,y,1)];
	for i = 2:n-1;
		yi = CDS(x,y,i);
		dy(end+1) = yi;
	end
	dy(end+1) = BDS(x,y,n);

	dy_real = df_list{fi}(x);
	e = rms(dy, dy_real);
	err(end+1) = e;
end
average_error = sum(err)/length(err);
disp('average rms error for n=10 points for each function')
err



% approximate derivatives using spectral differentiation for 10 points
n = 10;
x = linspace(-1,1,n)';
err = [];
for fi = 1:length(f_list);
	y = f_list{fi}(x);
	D = DerivMatrix(x,n-1);
	dy = D*y;
	dy_real = df_list{fi}(x);
	e = rms(dy, dy_real);
	err(end+1) = e;
end
average_error = sum(err)/length(err);
disp('average rms error for n=10 points for each function')
err


% central differencing for n = 3: 30 points
avg_err = [];
ns = 4:2:30;
for n = ns;
	x = linspace(-1,1,n);
	err = [];
	for fi = 1:length(f_list);
		y = f_list{fi}(x);
		dy = [FDS(x,y,1)];
		for i = 2:n-1;
			yi = CDS(x,y,i);
			dy(end+1) = yi;
		end
		dy(end+1) = BDS(x,y,n);

		dy_real = df_list{fi}(x);
		e = rms(dy, dy_real);
		err(end+1) = e;
	end
	average_error = sum(err)/length(err);
	avg_err(end+1) = average_error;

end

plot(ns,avg_err)
title('differencing derivative error for number of points n')
xlabel('n (number of divisions)')
ylabel('average rms error')


% central differencing for n = 3: 30 points
avg_err = [];
ns = 4:2:30;
for n = ns;
	x = linspace(-1,1,n)';
	err = [];
	for fi = 1:length(f_list);

		y = f_list{fi}(x);
		D = DerivMatrix(x,n-1);
		dy = D*y;
		dy_real = df_list{fi}(x);
		e = rms(dy, dy_real);
		err(end+1) = e;
	end
	average_error = sum(err)/length(err);
	avg_err(end+1) = average_error;

end

figure
plot(ns,avg_err)
title('spectral derivative error for number of points n')
xlabel('n (number of divisions)')
ylabel('average rms error')
disp('spectral error seems to decrease until n=14, then begins to rise')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Functions used for making numerical differentiation %%
function res = f1(x)
	res = abs(x.^3);
end

function res = f2(x)
	res = exp(-x.^(-2));
end

function res = f3(x)
	res = 1./(1+x.^2);
end

function res = f4(x)
	res = x.^10;
end

%%% exact derivative functions %%%
function res = df1(x)
	res = 3*abs(x.^2);
end

function res = df2(x)
	%res = 2*exp(-x.^(-2))./x.^3;
	res = 2*f2(x)./x.^3;
end

function res = df3(x)
	%res = -2.*x.*(1./(1+x.^2)).^2;
	res = -2.*x.*f3(x).^2;
end

function res = df4(x)
	res = 10*x.^9;
end

%% Functions for differencing schemes %%
function d = FDS(x,y,i)
	% Return 1st derivative approximation (forward differencing)
	delta = x(2) - x(1);
	num = y(i+1) - y(i);
	d = num/delta;
end


function d = BDS(x,y,i)
	% Return 1st derivative approximation (backward differencing)
	delta = x(2) - x(1);
	num = y(i) - y(i-1);
	d = num/delta;
end


function d = CDS(x,y,i)
	% Return 1st derivative approximation (central differencing)
	delta = x(2) - x(1);
	num = y(i) - y(i-1);
	d = num/delta;
end

%% other helping functions %%

function err = rms(a1, a2)
	% estimate error with root mean square between 2 arrays
	err = sqrt(sum((abs(a1 - a2)).^2)/length(a1));
end

function D = DerivMatrix(x,n)
	% Find D matrix (given in lecture)
	D = zeros(n+1,n+1);
	for i = 1:n+1;
		num(i)=1.0;
		for k=1:n+1
			if k~=i
				num(i) = num(i)*(x(i)-x(k));
			end
		end
	end
	for j=1:n+1;
		den(j)=1.0;
		for k =1:n+1
			if k~=j
				den(j) = den(j)*(x(j)-x(k));
			end
		end
	end
	for i=1:n+1
		for j=1:n+1
			if i~=j
				D(i,j)=num(i)/(den(j)*(x(i)-x(j)));
			end
		end
	end
	for i=1:n+1
		for k=1:n+1
			if i~=k
				D(i,i) = D(i,i) + 1./(x(i)-x(k));
			end
		end
	end
end
