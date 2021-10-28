%% Clear
clf
clear
format compact

%% Get images
% sequenceName = "JPEGS_FROM_LISAM";
% sequenceName = "MEMORIAL";
% Images below come from the same dataset and all have some weird
% horizontal lines artifacts that are primarily visible in darger regions
% sequenceName = "PARKING";
% sequenceName = "CONSTRUCTION";
sequenceName = "COUNTER";

[pixelValues, exposureTimes] = getImageSequence(sequenceName);
[R, G, B] = extractRGB(pixelValues);

disp(strcat("Number of images: ", num2str(length(pixelValues))));
% montage(pixelValues);

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

%% Tone-mapping
radMap = (cat(3, irradianceR, irradianceG, irradianceB));
trueToneMapped = tonemap(radMap);

%% Durand
contrast = 10; % Since wse want a contrast of 1:100 that LDR-images can handlec
luminanceMap = (radMap(:,:,1) * 20 + radMap(:,:,2) * 40 + radMap(:,:,3) * 1) / 61; % Intensity map
rScaled = radMap(:,:,1) ./ luminanceMap;
gScaled = radMap(:,:,2) ./ luminanceMap;
bScaled = radMap(:,:,3) ./ luminanceMap;
logLuminance = log10(luminanceMap);

baseLayer = imbilatfilt(logLuminance, 0.4, 0.02*size(logLuminance, 2)); % This is very slow for larger images

detailLayer = logLuminance - baseLayer; % High-pass filtered
maxBaseVal = max(baseLayer(:));
minBaseVal = min(baseLayer(:));

gamma = log10(contrast) / (maxBaseVal - minBaseVal); % To scale base-layer
logOutputMap = baseLayer * gamma + detailLayer;
outputMap = 10 .^ logOutputMap;
image(:,:,1) = outputMap .* rScaled; % Re-apply colors
image(:,:,2) = outputMap .* gScaled; % Re-apply colors
image(:,:,3) = outputMap .* bScaled; % Re-apply colors
scale = 1.0 / (10 ^ (maxBaseVal * gamma));
imageDurand = min(1.0, max(0.0, power(image*scale, 1.0/2.2)));

% subplot(1,2,1);
% imshow(image);
% title('Bilateral, Dorsey')
% subplot(1,2,2);
% imshow(trueToneMapped);
% title('MATLAB Tonemap')

%% Drago
% XYZ-transformation as defined in the paper
% Can't use MATLAB's rgb2xyz, probably related to the color space
% Rec. 709 used by Drago which is unavailable in rgb2xyz
xyz(:,:,1) = 0.412453 .* radMap(:,:,1) + 0.357580 .* radMap(:,:,2) + 0.180423 .* radMap(:,:,3);
xyz(:,:,2) = 0.212671 .* radMap(:,:,1) + 0.715160 .* radMap(:,:,2) + 0.072169 .* radMap(:,:,3);
xyz(:,:,3) = 0.019334 .* radMap(:,:,1) + 0.119193 .* radMap(:,:,2) + 0.950227 .* radMap(:,:,3);

% Convert to Yxy/xyY
W = sum(xyz,3);
Yxy(:,:,1) = xyz(:,:,2);        % Y
Yxy(:,:,2) = xyz(:,:,1) ./ W;	% x
Yxy(:,:,3) = xyz(:,:,2) ./ W;	% y

% Calculate some useful stats
numElements = numel(Yxy(:,:,1));
maxLuminance = max(max(Yxy(:,:,1)));
logLuminanceSum = sum(log(reshape(Yxy(:,:,1), [1 numElements] )));
logAvgLuminance = logLuminanceSum / numElements;
avgLuminance = exp(logAvgLuminance);
maxLuminanceW = (maxLuminance / avgLuminance); % Scaled max luminance

% This value has worked well
b = 0.85;

% Bias function
bT = (Yxy(:,:,1) ./ maxLuminanceW) .^ ( log(b) / log(0.5) );

% Calculate new luminance values
coeff = (100 * 0.01) / log10(maxLuminanceW + 1);
Yxy(:,:,1) = Yxy(:,:,1) ./ avgLuminance;
Yxy(:,:,1) = ( log(Yxy(:,:,1) + 1) ./ log(2 + bT .* 8) ) .* coeff;

% Convert Yxy/xyY to XYZ
newW = Yxy(:,:,1) ./ Yxy(:,:,3);
xyz(:,:,2) = Yxy(:,:,1);
xyz(:,:,1) = newW .* Yxy(:,:,2);
xyz(:,:,3) = newW -xyz(:,:,1) - xyz(:,:,2);

% XYZ to RGB
image(:,:,1) = 3.240479 .* xyz(:,:,1) + -1.537150 .* xyz(:,:,2) + -0.498535 .* xyz(:,:,3);
image(:,:,2) = -0.969256 .* xyz(:,:,1) + 1.875992 .* xyz(:,:,2) + 0.041556 .* xyz(:,:,3);
image(:,:,3) = 0.055648 .* xyz(:,:,1) + -0.204043 .* xyz(:,:,2) + 1.057311 .* xyz(:,:,3);

% correct gamma
imageDrago = fixGamma(image, 2.7);

%% global tonemapping
%osäker på om det behövs samplas igen
% subplot(1,1,1);
hsv = rgb2hsv(radMap); %transforms the image values
luminance=luminanceMap;
chrominance=hsv(:,:,2); %chrominance = saturation channel
[rows, cols]=size(radMap(:,:,1));
T=rows*cols; %total amount of pixels
logLum=log10(luminance);
deltaB=(max(max(logLum)) - min(min(logLum)) )/100;

displayMin=100; % ?
displayMax=100000000000000000000000000; % ?
logdmin=log10(displayMin);
logdmax=log10(displayMax);
bincount=100;
h=imhist(luminance,bincount) / T;
logDisplay=zeros(size(luminance));

%cumsum = P(b)
%for k=1:5:rows
   % for m=1:5:cols
      %  logDisplay(k,m)=logdmin + (logdmax - logdmin).*cumsum(logLum(k,m)/T);
      %  lumDis=logDisplay./luminance;
      %  %lumDisDerivative=diff(logDisplay, luminance);
   % end
%end

fb=imhist(logLum, bincount);

for i=1:1:bincount
    if fb(i) > (T*deltaB) / (logdmax - logdmin)
            fb(i) = (T*deltaB) / (logdmax - logdmin);
    end
    
    if fb(i) < displayMin
            fb(i)=T;
    end
end

globalToneMappedImage=histeq(luminance,fb);
%globalToneMappedImage=luminance.*fb;
%imshow(luminance);
%imshow(lumDis);
%imshow(logDisplay);
%imshow(globalToneMappedImage);

imageGlobal(:,:,1) = globalToneMappedImage .* rScaled.*0.7; % blev för rött så ändrade skalningen
imageGlobal(:,:,2) = globalToneMappedImage .* gScaled; % Re-apply colors
imageGlobal(:,:,3) = globalToneMappedImage .* bScaled; % Re-apply colors
imshow(imageGlobal);

%imhist(luminance,bincount)
imhist(globalToneMappedImage,bincount)
%imhist(logLum, bincount);
%imhist(logDisplay,bincount);

%% Another global one
a = 0.58;
saturation = 0.7;
imshow(reinhardGlobal(irradianceR, irradianceG, irradianceB, a, saturation));

%% And a third global one (very simple)
gamma = 0.2;
RGB = gammaToneMapping(irradianceR, irradianceG, irradianceB, gamma);
imshow(RGB)

%% Adaptation of the gamma one to luminance
gamma = 0.2;
saturation = 0.4;
RGB = gammaLuminanceToneMapping(irradianceR, irradianceG, irradianceB, ...
    gamma, saturation);
imshow(RGB)

%%

subplot(2,2,1);
imshow(trueToneMapped);
title('MATLAB')
subplot(2,2,2);
imshow(imageDurand);
title('Durand')
subplot(2,2,3);
imshow(imageDrago);
title('Drago')
subplot(2,2,4);
imshow(imageGlobal);
title('Global')