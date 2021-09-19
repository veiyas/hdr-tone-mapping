function [pixelValues, exposureTimes] = getImageSequence(sequenceId)

dataDirectory = "../data/";
pixelValues = {};
exposureTimes = {};

if strcmp(sequenceId, "JPEGS_FROM_LISAM")
    ids = 5165:5178;
    for i = 1:length(ids)
        name = strcat(dataDirectory, "/HDR/IMG_", num2str(ids(i)), ".JPG");
        pixelValues{i} = imread(name);
        exposureTimes{i} = imfinfo(name).DigitalCamera.ExposureTime;
    end
elseif strcmp(sequenceId, "MEMORIAL")
    files = [...
        "memorial0061.png" 0.03125
        "memorial0062.png" 0.0625
        "memorial0063.png" 0.125
        "memorial0064.png" 0.25
        "memorial0065.png" 0.5
        "memorial0066.png" 1
        "memorial0067.png" 2
        "memorial0068.png" 4
        "memorial0069.png" 8
        "memorial0070.png" 16
        "memorial0071.png" 32
        "memorial0072.png" 64
        "memorial0073.png" 128
        "memorial0074.png" 256
        "memorial0075.png" 512
        "memorial0076.png" 1024];

    for i = 1:size(files, 1)
        name = strcat(dataDirectory, "memorial/", files(i,1));
        image = imread(name);
        % Crop the stupid blue borders
        pixelValues{i} = image(48:(768-24), 4:(512-16), :);
        exposureTimes{i} = 1 / str2double(files(i,2));
    end
end

end

