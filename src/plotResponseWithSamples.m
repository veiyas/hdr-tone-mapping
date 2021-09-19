function plotResponseWithSamples(lE, logExp, Z, g)
markerColor = [1 0.5 1];
scatter(lE + logExp, Z, 'x', 'MarkerEdgeColor', markerColor,...
    'LineWidth', 0.1)
hold on
xlabel("log exposure")
ylabel("pixel value")
plot(g, 0:255, 'k', 'LineWidth', 2)
hold off
end

