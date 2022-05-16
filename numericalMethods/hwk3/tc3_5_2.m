%example from workshop 8.7

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




function D = DerivMatrix(x,n)
	% Spectral differentiation
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
