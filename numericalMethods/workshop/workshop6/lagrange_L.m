function L = lagrange_L(i,nodes,x)
	% i - ith lagrange polynomial to be calculated
	% nodes: known x vals of the original function
	% x_range domain (many points for interpolation)

L = ones(length(x),1);

for k = 1:length(x);
    for j = 1:length(nodes);
        if j~=i
            L(k) = L(k)*(x(k)-nodes(j))/(nodes(i)-nodes(j));
        end
    end
end

end
