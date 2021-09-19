%% Clear
clf
clear
format compact

%% Get images
sequenceName = "JPEGS_FROM_LISAM";
% sequenceName = "MEMORIAL";

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

lambda = 150;
% This feels kinda weird, comments in gSolve for B dont match code?
logExp = repmat(logExposureTimes, numPixelSamples)';

[gRed,~]   = gSolve(ZRed,   logExp, lambda, w);
[gGreen,~] = gSolve(ZGreen, logExp, lambda, w);
[gBlue,~]  = gSolve(ZBlue,  logExp, lambda, w);

%% Plot the recovered camera response function
plot(gRed, 0:255, 'r')
xlabel("log exposure")
ylabel("pixel value")
grid on
hold on
plot(gGreen, 0:255, 'g')
plot(gBlue, 0:255, 'b')
hold off

%% Constructing the radiance map
% A different weight function can be used... see (10.10) in computer vision

[logIrradianceR] = constructRadianceMap(R, w, logExposureTimes, gRed);
[logIrradianceG] = constructRadianceMap(G, w, logExposureTimes, gGreen);
[logIrradianceB] = constructRadianceMap(B, w, logExposureTimes, gBlue);

%% Visualizing the radiance maps with color scale
colormap jet
subplot(221); imagesc(logIrradianceR); colorbar; title("R")
subplot(222); imagesc(logIrradianceG); colorbar; title("G")
subplot(223); imagesc(logIrradianceB); colorbar; title("B")
