function [score] = score_image(I)
    score = 0;
    
    R = double(I(:,:,1));
    G = double(I(:,:,2));
    B = double(I(:,:,3));
    
    S = size(R);
    N = S(1);
    M = S(2);
    
    
    for x = 1:M
        for y = 1:N
            score = score + abs(R(y,x) - G(y,x))^2 + abs(B(y,x) - G(y,x))^2 + abs(R(y,x) - B(y,x))^2;
        end
    end
end

