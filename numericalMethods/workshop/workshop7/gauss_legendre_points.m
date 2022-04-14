function [w,x] = gauss_legendre_points(k)
    
    syms t
    x = double(vpasolve(legendreP(k,t) == 0));
    
    w = 2*(1-x.^2)./((k+1)^2.*(legendreP(k+1,x).^2));
    
end