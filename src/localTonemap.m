function [result] = localTonemap(logLum, baseLayer, detailScale, plotLayers)
detailLayer = logLum - baseLayer;
detailLayerScaled = detailScale .* detailLayer;
result = exp(baseLayer + detailLayerScaled);

if plotLayers
    % Plot filtered layers and result
    figure;
    subplot(2,2,1)
    imshow(logLum);
    subplot(2,2,2)
    imshow(baseLayer);
    subplot(2,2,3);
    imshow(detailLayer);
    subplot(2,2,4);
    imshow(result); 
end
end

