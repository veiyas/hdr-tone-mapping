function [R, G, B] = extractRGB(images)
% Assumes same size for all images
width = size(images{1}, 1);
height = size(images{1}, 2);

R = zeros(width, height, length(images));
G = zeros(width, height, length(images));
B = zeros(width, height, length(images));

for j=1:length(images)
   R(:,:,j) = images{j}(:,:,1);
   G(:,:,j) = images{j}(:,:,2);
   B(:,:,j) = images{j}(:,:,3);
end
end

