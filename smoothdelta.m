function [y] = smoothdelta(x)
    x = min(1, max(0, x));
    
    y = 6.*x - 6.*x.^2;
end

