function [pixelValues, exposureTimes] = getImageSequence(sequenceId)

dataDirectory = "../data/";

if strcmp(sequenceId, "JPEGS_FROM_LISAM")
    ids = 5165:5178;
    pixelValues = {};
    exposureTimes = {};
    for i = 1:length(ids)
        name = strcat(dataDirectory, "/HDR/IMG_", num2str(ids(i)), ".JPG");
        pixelValues{i} = imread(name);
        exposureTimes{i} = imfinfo(name).DigitalCamera.ExposureTime;
    end
end 

end

