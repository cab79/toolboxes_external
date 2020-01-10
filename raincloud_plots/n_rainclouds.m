%% simultaneous plot of N rainclouds - 
% plots density plot, boxplot, and 1-d scatter concurrently for several
% data vectors
% Use as h = n_rainclouds(X, cl), 
% where X is an array of N cells and cl is an N x 3 array of RGB values. 
% h is a cell array of handles for the various figure parts.
% Written by Tom Marshall. www.tomrmarshall.com

function h = n_rainclouds(X, cl)

% clouds
for i = 1:length(X)
    
    [f{i}, xi{i}] = ksdensity(X{i});
    % density plot
    h.dx{i} = area(xi{i}, f{i}); hold on
    set(h.dx{i}, 'FaceColor', cl(i,:));
    set(h.dx{i}, 'EdgeColor', 'none');
    set(h.dx{i}, 'FaceAlpha', 0.5);
    
end

% make some space under the clouds for the drops
yl = get(gca, 'YLim');
set(gca, 'YLim', [-yl(2) yl(2)]);

% width of boxplot / drops
wdth = yl(2)*0.5;

% raindrops
for i = 1:length(X)
    
    % jitter for raindrops
    jit{i} = (rand(size(X{i})) - 0.5) * wdth;
    
    % raindrops
    h.dr{i} = scatter(X{i}, jit{i} - yl(2)/2);
    h.dr{i}.SizeData = 32;
    h.dr{i}.MarkerFaceColor = cl(i,:);
    h.dr{i}.MarkerEdgeColor = 'none';
    set(h.dr{i},'MarkerFaceAlpha', 0.5);
    
end

% mean and quantiles
for i = 1:length(X)
    
    Y{i} = quantile(X{i}, [0.25 0.75 0.5 0.02 0.98]);
    offset(i) = (wdth * 0.1) * (i-1) * length(X);
    h.bx{i,1} = line([Y{i}(1) Y{i}(2)], wdth - [offset(i) offset(i)], 'col', 'k', 'LineWidth', 2);
    h.bx{i,2} = scatter(mean(X{i}),wdth - offset(i), 120, cl(i,:));
    h.bx{i,2}.MarkerFaceColor = cl(i,:);
    h.bx{i,2}.MarkerEdgeColor = [0 0 0];
    h.bx{i,2}.LineWidth = 2;
    
end
