function Ni = newton_N(x,i,xpoints)

Ni = ones(size(x));

for l=1:i
    Ni = Ni.*(x-xpoints(l));
end

end