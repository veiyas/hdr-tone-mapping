%% Clear
clf
clear
format compact

%% Get images
[pixelValues, exposureTimes] = getImageSequence("JPEGS_FROM_LISAM");
[R, G, B] = extractRGB(pixelValues);
disp(strcat("Number of images: ", num2str(length(pixelValues))));
montage(pixelValues);

%% Not sure if this is correct yet
logExposureTimes = log(cell2mat(exposureTimes))';

% Since using all pixels would produce collossal matrices later, the
% images has to be sampled in some way.
numPixelSamples = 150;
ZRed   = zeros(numPixelSamples, length(logExposureTimes));
ZGreen = zeros(numPixelSamples, length(logExposureTimes));
ZBlue  = zeros(numPixelSamples, length(logExposureTimes));

for j=1:length(pixelValues)
    tempR = R(:,:,j);
    tempG = G(:,:,j);
    tempB = B(:,:,j);
    % Just choosing random samples is probably not optimal
    % So thats something we could improve
    ZRed(:,j)   = datasample(tempR(:), numPixelSamples);
    ZGreen(:,j) = datasample(tempG(:), numPixelSamples);
    ZBlue(:,j)  = datasample(tempB(:), numPixelSamples);
end

% NOTE Maybe this creates off by 1 when used in gSolve and/or
% constructRadianceMap, but also maybe not; will investigate
w = arrayfun(@weightingFunction, 0:255);

lambda = 500; % Not sure what values are reasonable for this
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
subplot(221); imagesc(logIrradianceR, [-1.5, 7]); colorbar; title("R")
subplot(222); imagesc(logIrradianceG, [-1.5, 7]); colorbar; title("G")
subplot(223); imagesc(logIrradianceB, [-1.5, 7]); colorbar; title("B")
