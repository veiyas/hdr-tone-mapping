%% Clear
clf
clear
format compact

%% Get images
% sequenceName = "JPEGS_FROM_LISAM";
sequenceName = "MEMORIAL";

[pixelValues, exposureTimes] = getImageSequence(sequenceName);
[R, G, B] = extractRGB(pixelValues);
% [R, G, B] = extractLAB(pixelValues);

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
% plot(gRed, 0:255, 'r')
% xlabel("log exposure")
% ylabel("pixel value")
% grid on
% hold on
% plot(gGreen, 0:255, 'g')
% plot(gBlue, 0:255, 'b')
% hold off

%% Plot camera response with the samples
% subplot(221); plotResponseWithSamples(lERed, logExp, ZRed, gRed);
%     title('R');
% subplot(222); plotResponseWithSamples(lEGreen, logExp, ZGreen, gGreen);
%     title('G');
% subplot(223); plotResponseWithSamples(lEBlue, logExp, ZBlue, gBlue);
%     title('B');

%% Constructing the radiance map
% A different weight function can be used... see (10.10) in computer vision

[irradanceR] = constructRadianceMap(R, w, logExposureTimes, gRed);
[irradianceG] = constructRadianceMap(G, w, logExposureTimes, gGreen);
[irradianceB] = constructRadianceMap(B, w, logExposureTimes, gBlue);

%% Visualizing the radiance maps with color scale
colormap jet
subplot(221); imagesc(irradanceR); colorbar; title("R")
subplot(222); imagesc(irradianceG); colorbar; title("G")
subplot(223); imagesc(irradianceB); colorbar; title("B")

%% Construct CIE-Lab image
rgbIm = (cat(3, irradanceR, irradianceG, irradianceB));
% rgbIm = im2double(rgbIm);

labIm = rgb2lab(rgbIm);
logLum = log(labIm(:,:,1));

% rbgImMin = min(rgbIm(:))
% rbgImMax = max(rgbIm(:))

%% Local filtering
plotLayers = 0; % True or false
baseLayerGauss = imgaussfilt(logLum, 10); % param2: Standard deviation
baseLayerBilat = imbilatfilt(logLum, 100); % param2: degrees of smoothing

luminanceScale = 0.8;
resultGauss = localTonemap(logLum, baseLayerGauss, luminanceScale, plotLayers);
resultBilat = localTonemap(logLum, baseLayerBilat, luminanceScale, plotLayers);

rgbGauss = lab2rgb(cat(3, resultGauss, labIm(:,:,2), labIm(:,:,3)));
rgbBilat = lab2rgb(cat(3, resultBilat, labIm(:,:,2), labIm(:,:,3)));

% Plot tonemapped and non-tonemapped
if plotLayers
   figure;
end
subplot(1,3,1);
imshow(rgbIm);
title('No tonemapping')
subplot(1,3,2);
imshow(rgbGauss);
title('Gauss low-pass')
subplot(1,3,3);
imshow(rgbBilat);
title('Bilateral')





