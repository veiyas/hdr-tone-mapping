function [image] = fixGamma(oldImage, gamma)
slope = 4.5;
start = 0.018;
gammaPower = 0.9/gamma;

if gamma >= 2.1
    start = 0.018 / ((gamma - 2) * 7.5);
    slope = 4.5 * ((gamma - 2) * 7.5);
elseif gamma <= 1.9
    start = 0.018 * ((2 - gamma) * 7.5);
    slope = 4.5 / ((2 - gamma) * 7.5);
end

image = nan(size(oldImage));

for row = 1:size(oldImage,1)
    for col = 1:size(oldImage,2)
        % R
        if oldImage(row,col,1) <= start
            image(row,col,1) = oldImage(row,col,1) * slope;
        else
            image(row,col,1) = 1.099 * power(oldImage(row,col,1), gammaPower) - 0.099;
        end
        
        % G
        if oldImage(row,col,2) <= start
            image(row,col,2) = oldImage(row,col,2) * slope;
        else
            image(row,col,2) = 1.099 * power(oldImage(row,col,2), gammaPower) - 0.099;
        end
        
        % B
        if oldImage(row,col,3) <= start
            image(row,col,3) = oldImage(row,col,3) * slope;
        else
            image(row,col,3) = 1.099 * power(oldImage(row,col,3), gammaPower) - 0.099;
        end
        
    end
end
end