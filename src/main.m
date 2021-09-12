%% Clear
clear
format compact

%% Get images
[pixelValues, exposureTimes] = getImageSequence("JPEGS_FROM_LISAM");
disp(strcat("Number of images: ", num2str(length(pixelValues))));
montage(pixelValues);
