%% raincloud_plot - plots a combination of half-violin, boxplot,  and raw
% datapoints (1d scatter).
% Use as h = raincloud_plot(X, cl), where X is a data vector and cl is an
% RGB value. h is a cell array of handles for the various figure parts.
% Optional 3rd input argument 'density_type' can be 'ks' (default) or 'rash'.
% 'ks' uses matlab's 'ksdensity'. 'rash' uses 'rst_RASH' from Cyril
% Pernet's robust statistics toolbox (which must be on the matlab path).
% Based on https://micahallen.org/2018/03/15/introducing-raincloud-plots/
% Inspired by https://m.xkcd.com/1967/
% Written by Tom Marshall. www.tomrmarshall.com
% Thanks to Jacob Bellmund for some improvements


function h = raincloud_plot(X, cl, density_type, box_on)

if ~exist('density_type', 'var') | isempty(density_type)
    density_type = 'ks'; % default is 'ks', can also be 'rash'
end

if nargin < 4
    box_on = 1;
end

% calculate kernel density
switch density_type
    case 'ks'
        [f, Xi] = ksdensity(X);
    case 'rash'
        try
            [Xi, f] = rst_RASH(X);
        catch
            disp('you''ve specified density_type = ''RASH'', but something''s gone wrong.')
            disp('Have you downloaded Cyril Pernet''s robust stats toolbox?');
        end
end

% density plot
h{1} = area(Xi, f); hold on
set(h{1}, 'FaceColor', cl);
set(h{1}, 'EdgeColor', [0.1 0.1 0.1]);
set(h{1}, 'LineWidth', 2);

% make some space under the density plot for the boxplot
yl = get(gca, 'YLim');
set(gca, 'YLim', [-yl(2) yl(2)]);

% width of boxplot
wdth = yl(2)*0.5;

% jitter for raindrops
jit = (rand(size(X)) - 0.5) * wdth;

% info for making boxplot
Y = quantile(X, [0.25 0.75 0.5 0.02 0.98]);

% raindrops
h{2} = scatter(X, jit - yl(2)/2);
h{2}.SizeData = 10;
h{2}.MarkerFaceColor = cl;
h{2}.MarkerEdgeColor = 'none';

if box_on
    % 'box' of 'boxplot'
    h{3} = rectangle('Position', [Y(1) -yl(2)/2-(wdth*0.5) Y(2)-Y(1) wdth]);
    set(h{3}, 'EdgeColor', 'k')
    set(h{3}, 'LineWidth', 2);
    % could also set 'FaceColor' here as Micah does, but I prefer without
    
    % mean line
    h{4} = line([Y(3) Y(3)], [-yl(2)/2-(wdth*0.5) -yl(2)/2+(wdth*0.5)], 'col', 'k', 'LineWidth', 2);
    
    % whiskers
    h{5} = line([Y(2) Y(5)], [-yl(2)/2 -yl(2)/2], 'col', 'k', 'LineWidth', 2);
    h{6} = line([Y(1) Y(4)], [-yl(2)/2 -yl(2)/2], 'col', 'k', 'LineWidth', 2);
end
