hdr(:,:,1) = irradianceR;
hdr(:,:,2) = irradianceG;
hdr(:,:,3) = irradianceB;

Emax = max([irradianceR; irradianceG; irradianceB], [], 'all');
Emin = min([irradianceR; irradianceG; irradianceB], [], 'all');

decidedMin = Emin;
decidedMax = Emax;

Enorm(:,:,1) = (irradianceR - decidedMin) / (decidedMax - decidedMin);
Enorm(:,:,2) = (irradianceG - decidedMin) / (decidedMax - decidedMin);
Enorm(:,:,3) = (irradianceB - decidedMin) / (decidedMax - decidedMin);

minn = 0.0;
maxx = 0.001;

lin = (Enorm - minn) / (maxx - minn);
Enorm = max(min(lin, 1), 0);

imshow(Enorm)

%%

% log histogram:
% https://se.mathworks.com/matlabcentral/answers/102406-how-can-i-plot-a-histogram-with-a-logarithmic-x-axis
[~,edges] = histcounts(log10(hdr));
histogram(hdr,10.^edges)
set(gca, 'xscale','log')
axis tight

%%
[~,edges] = histcounts(log10(Enorm));
histogram(Enorm,10.^edges)
set(gca, 'xscale','log')
axis tight

%%
x = linspace(0,1);
plot(x, x, 'LineWidth', 2);
xlabel('Irradiance, E');
ylabel('Pixel intensity');
xticks([0 1])
xticklabels({'Emin','Emax'})
yticks([0 1])
axis equal
axis tight

%%
x = linspace(0,1,1000);
minn = 0.2;
maxx = 0.5;
line = (x - minn) ./ (maxx - minn);
y = max(min(line, 1), 0);
plot(x, y, 'LineWidth', 2);
xlabel('Irradiance, E');
ylabel('Pixel intensity');
xticks([0 1])
xticklabels({'Emin','Emax'})
yticks([0 1])
axis equal
axis tight