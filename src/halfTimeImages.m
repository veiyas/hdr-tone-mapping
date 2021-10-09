%% Clear
clf
clear
format compact

%sequenceName = "JPEGS_FROM_LISAM";
sequenceName = "MEMORIAL";

[pixelValues, exposureTimes] = getImageSequence(sequenceName);
[R, G, B] = extractRGB(pixelValues);
disp(strcat("Number of images: ", num2str(length(pixelValues))));
montage(pixelValues, 'Size', [2 7]);

%% Not sure if this is correct yet
logExposureTimes = log(cell2mat(exposureTimes))';

% Since using all pixels would produce collossal matrices later, the
% images has to be sampled in some way.
numPixels = size(pixelValues{1}, 1) * size(pixelValues{1}, 2);
numPixelSamples = 200;
ZRed   = zeros(numPixelSamples, length(logExposureTimes));

step = numPixels / numPixelSamples;
sampleIndices = floor((1:step:numPixels));

for j=1:length(exposureTimes)
    tempR = reshape(R(:,:,j), numPixels, 1);
    ZRed(:,j)   = tempR(sampleIndices);
end

w = arrayfun(@weightingFunction, 0:255);
%w = ones(256, 1);

lambda = 150;
% This feels kinda weird, comments in gSolve for B dont match code?
logExp = repmat(logExposureTimes, 1, numPixelSamples)';

[gRed, lERed] = gSolve(ZRed, logExp, lambda, w);

plotResponseWithSamples(lERed, logExp, ZRed, gRed);

%%
xlabel("log exposure")
ylabel("pixel value")
plot(gRed, 0:255, 'k', 'LineWidth', 2)

%%
[logIrradianceR] = constructRadianceMap(R, w, logExposureTimes, gRed);
imagesc(logIrradianceR); axis equal; colorbar

%%
plot(0:255, w);
xlabel("pixel value")
ylabel("w")