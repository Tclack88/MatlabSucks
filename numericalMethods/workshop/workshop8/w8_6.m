% derivative approximations using central differencing schemes
% Applied to thermoclines in a lake

z = [ 0 2.3 4.9 9.1 13.7 18.3 22.9 27.2]; % depth of lake (m)
T = [22.8 22.8 22.8 20.6 13.9 11.7 11.1 11.1]; % T(C)

% Thermoclines occur at points where the 2nd derivative is 0
% Since the points of z are not evenly spaced, we must use a generic calculation
% to find the ith 2nd derivative

d2z = zeros(1,length(z));

d2z(1) = FDS_2df(z,T,1);
d2z(end) = BDS_2df(z,T,length(z));
for i = 2:length(z)-1;
	d2z(i) = CDS_2df(z,T,i);
end

plot(d2z, z)
hold on
xline(0,'--k')
ylabel('depth')
xlabel('2nd derivative = 0')
title('Thermoclines in a lake')


%%%%%%%%%%%%% functions %%%%%%%%%%%%%%%%%

function d2 = CDS_2df(x,y , i)
	num = y(i+1)*(x(i)-x(i-1)) - y(i)*(x(i+1)-x(i-1)) + y(i-1)*(x(i+1)-x(i));
	denom = (x(i+1)-x(i-1))*(x(i+1)-x(i))*(x(i)-x(i-1))/2;
	d2 = num/denom;
end


function d2 = FDS_2df(x,y,i)
	delta = x(i+1) - x(i);
	num = y(i+2) - 2*y(i+1) + y(i);
	d2 = num/delta^2;
end


function d2 = BDS_2df(x,y,i)
	delta = x(i) - x(i-1);
	num = y(i) - 2*y(i-1) + y(i-2);
	d2 = num/delta^2;
end
