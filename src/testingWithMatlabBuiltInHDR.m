ids = 5165:5178;
fileNames = {};
exposureTimes = {};
for i = 1:length(ids)
    fileNames{i} = strcat("../data//HDR/IMG_", num2str(ids(i)), ".JPG");
    exposureTimes{i} = imfinfo(fileNames{i}).DigitalCamera.ExposureTime;
end

%%
crf = camresponse(fileNames);

range = 0:length(crf)-1;

hold on
plot(crf(:,1),range,'--r','LineWidth',2);
plot(crf(:,2),range,'-.g','LineWidth',2);
plot(crf(:,3),range,'-.b','LineWidth',2);
xlabel('Log-Exposure');
ylabel('Image Intensity');
title('Camera Response Function');
grid on
axis('tight')
legend('R-component','G-component','B-component','Location','southeast')
hold off

%%
expTimes = cell2mat(exposureTimes);
hdr = makehdr(fileNames, 'RelativeExposure', expTimes./expTimes(1));

%%
imagesc(hdr); colorbar

%%
rgb = tonemap(hdr);
imshow(rgb)