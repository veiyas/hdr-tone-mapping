%% Clear
clf
clear
format compact

%% Get images
[pixelValues, exposureTimes] = getImageSequence("JPEGS_FROM_LISAM");
disp(strcat("Number of images: ", num2str(length(pixelValues))));
montage(pixelValues);

%% Not sure if this is correct yet
logExposureTimes = log(cell2mat(exposureTimes))';

% Since using all pixels would produce collossal matrices later, the
% images has to be sampled in some way.
% numPixels = size(pixelValues{1}, 1) * size(pixelValues{1}, 2);
numPixelSamples = 150;
ZRed   = zeros(numPixelSamples, length(logExposureTimes));
ZGreen = zeros(numPixelSamples, length(logExposureTimes));
ZBlue  = zeros(numPixelSamples, length(logExposureTimes));

for i=1:length(pixelValues)
    R = pixelValues{i}(:,:,1);
    G = pixelValues{i}(:,:,2);
    B = pixelValues{i}(:,:,3);
    
    % Just choosing random samples is probably not optimal
    % So thats something we could improve
    ZRed(:,i)   = datasample(R(:), numPixelSamples);
    ZGreen(:,i) = datasample(G(:), numPixelSamples);
    ZBlue(:,i)  = datasample(B(:), numPixelSamples);
end

w = arrayfun(@weightingFunction, 0:255);
lambda = 150; % Not sure what values are reasonable for this
% This feels kinda weird, comments in gSolve for B dont match code?
B = repmat(logExposureTimes, numPixelSamples)';

[gRed,logIrradianceRed] = gSolve(ZRed, B, lambda, w);
[gGreen,logIrradianceGreen] = gSolve(ZGreen, B, lambda, w);
[gBlue,logIrradianceBlue] = gSolve(ZBlue, B, lambda, w);

%% Plot the recovered camera response function
plot(gRed, 0:255, 'r')
xlabel("log exposure")
ylabel("pixel value")
hold on
plot(gGreen, 0:255, 'g')
plot(gBlue, 0:255, 'b')
hold off




