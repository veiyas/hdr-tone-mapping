%% Clear
clf
clear
format compact

%% Get images
% sequenceName = "JPEGS_FROM_LISAM";
sequenceName = "MEMORIAL";
% Images below come from the same dataset and all have some weird
% horizontal lines artifacts that are primarily visible in darger regions
% sequenceName = "PARKING";
% sequenceName = "CONSTRUCTION";
% sequenceName = "COUNTER";

[pixelValues, exposureTimes] = getImageSequence(sequenceName);
[R, G, B] = extractRGB(pixelValues);

disp(strcat("Number of images: ", num2str(length(pixelValues))));
montage(pixelValues);

%% Not sure if this is correct yet
logExposureTimes = log(cell2mat(exposureTimes))';

% Since using all pixels would produce collossal matrices later, the
% images has to be sampled in some way.
numPixels = size(pixelValues{1}, 1) * size(pixelValues{1}, 2);
numPixelSamples = 500;
ZRed   = zeros(numPixelSamples, length(logExposureTimes));
ZGreen = zeros(numPixelSamples, length(logExposureTimes));
ZBlue  = zeros(numPixelSamples, length(logExposureTimes));

step = numPixels / numPixelSamples;
sampleIndices = floor((1:step:numPixels));

for j=1:length(exposureTimes)
    tempR = reshape(R(:,:,j), numPixels, 1);
    tempG = reshape(G(:,:,j), numPixels, 1);
    tempB = reshape(B(:,:,j), numPixels, 1);
    
    ZRed(:,j)   = tempR(sampleIndices);
    ZGreen(:,j) = tempG(sampleIndices);
    ZBlue(:,j)  = tempB(sampleIndices);
end

% NOTE Maybe this creates off by 1 when used in gSolve and/or
% constructRadianceMap, but also maybe not; will investigate
w = arrayfun(@weightingFunction, 0:255);

lambda = 100;
% This feels kinda weird, comments in gSolve for B dont match code?
logExp = repmat(logExposureTimes, 1, numPixelSamples)';

[gRed, lERed] = gSolve(ZRed, logExp, lambda, w);
[gGreen, lEGreen] = gSolve(ZGreen, logExp, lambda, w);
[gBlue, lEBlue] = gSolve(ZBlue, logExp, lambda, w);

%% Plot the recovered camera response function
plot(gRed, 0:255, 'r')
xlabel("log exposure")
ylabel("pixel value")
grid on
hold on
plot(gGreen, 0:255, 'g')
plot(gBlue, 0:255, 'b')
hold off

%% Plot camera response with the samples
% subplot(221); plotResponseWithSamples(lERed, logExp, ZRed, gRed);
%     title('R');
% subplot(222); plotResponseWithSamples(lEGreen, logExp, ZGreen, gGreen);
%     title('G');
% subplot(223); plotResponseWithSamples(lEBlue, logExp, ZBlue, gBlue);
%     title('B');

%% Constructing the radiance map
% A different weight function can be used... see (10.10) in computer vision

[irradianceR] = constructRadianceMap(R, w, logExposureTimes, gRed);
[irradianceG] = constructRadianceMap(G, w, logExposureTimes, gGreen);
[irradianceB] = constructRadianceMap(B, w, logExposureTimes, gBlue);

%% Visualizing the radiance maps with color scale
colormap jet
subplot(221); imagesc(irradianceR); colorbar; title("R")
subplot(222); imagesc(irradianceG); colorbar; title("G")
subplot(223); imagesc(irradianceB); colorbar; title("B")

%% Tonemap, Durand & Dorsey
% Potential luminance map
% luminanceMap = log10(0.2125 * hdrIm(:,:,1) + 0.7152 * hdrIm(:,:,2) + 0.0722 * hdrIm(:,:,3));

radmap = (cat(3, irradianceR, irradianceG, irradianceB));
% radmap = im2double(imadjust(radmap, [0 1], []));
trueToneMapped = tonemap(radmap);

contrast = 10; % Since wse want a contrast of 1:100 that LDR-images can handlec
intensityMap = (radmap(:,:,1) * 20 + radmap(:,:,2) * 40 + radmap(:,:,3) * 1) / 61; % Intensity map
rScaled = radmap(:,:,1) ./ intensityMap;
gScaled = radmap(:,:,2) ./ intensityMap;
bScaled = radmap(:,:,3) ./ intensityMap;
logIntensity = log10(intensityMap);

baseLayer = imbilatfilt(logIntensity, 0.4, 0.02*size(logIntensity, 2)); % This is too slow for larger images

detailLayer = logIntensity - baseLayer; % High-pass filtered
maxBaseVal = max(baseLayer(:));
minBaseVal = min(baseLayer(:));

gamma = log10(contrast) / (maxBaseVal - minBaseVal); % To scale base-layer
%Offset = max(max(lBase)) * Scale;
%lOmap = lBase * Scale + lDetail - Offset;
logOutputMap = baseLayer * gamma + detailLayer;
outputMap = 10 .^ logOutputMap;
image(:,:,1) = outputMap .* rScaled; % Re-apply colors
image(:,:,2) = outputMap .* gScaled; % Re-apply colors
image(:,:,3) = outputMap .* bScaled; % Re-apply colors
%maxVal = max(max(max(image)));
%image = image / maxVal;
scale = 1.0 / (10 ^ (maxBaseVal * gamma));
image = min(1.0, max(0.0, power(image*scale, 1.0/2.2)));

subplot(1,2,1);
imshow(image);
title('Bilateral, Dorsey')
subplot(1,2,2);
imshow(trueToneMapped);
title('MATLAB Tonemap')
