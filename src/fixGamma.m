function [image] = fixGamma(oldImage, gamma)
slope = 4.5;
start = 0.018;
fgamma = (0.45/gamma)*2;

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
        %red
        if oldImage(row,col,1) <= start
            image(row,col,1) = oldImage(row,col,1) * slope;
        else
            image(row,col,1) = 1.099 * power(oldImage(row,col,1), fgamma) - 0.099;
        end
        
        %green
        if oldImage(row,col,2) <= start
            image(row,col,2) = oldImage(row,col,2) * slope;
        else
            image(row,col,2) = 1.099 * power(oldImage(row,col,2), fgamma) - 0.099;
        end
        
        %blue
        if oldImage(row,col,3) <= start
            image(row,col,3) = oldImage(row,col,3) * slope;
        else
            image(row,col,3) = 1.099 * power(oldImage(row,col,3), fgamma) - 0.099;
        end
        
    end
end
end