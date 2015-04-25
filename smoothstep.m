function [y] = smoothstep(x)
    x = min(1, max(0, x));
    
    y = x.^2.*(3 - 2.*x);
end

