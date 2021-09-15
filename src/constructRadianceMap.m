function [logIrradiance] = constructRadianceMap(images, w, logDt, g)
% Takes one color channel only
% images is the images with the third dimension denoting which image it is

width = size(images, 1);
height = size(images, 2);
numImages = size(images, 3);
numPixels = width * height;

% Reshape to be more like the notation in the paper
Z = reshape(images, [numPixels, numImages]);

num = zeros(numPixels, 1);
den = zeros(numPixels, 1);

% It is very possible that this could be done in a smart matlab way
% but this will do (at least for now)
for i=1:numPixels   
    for j=1:numImages
        num(i) = num(i) + (w(Z(i,j) + 1) .* (g(Z(i,j) + 1) - logDt(j)));
        den(i) = den(i) + w(Z(i,j) + 1);
    end
end

lnE = num ./ den;
logIrradiance = reshape(lnE, [width, height]);

end

