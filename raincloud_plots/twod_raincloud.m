%% two-d raincloud plot -
% plots kernel densities, 2-d scatterplot of data and prediction intervals
% use like:
% h = twod_raincloud(X, cl)
% where X is a vector of n cells with the data values (2 rows)
% and cl is a 3 x n matrix of RGB values
% 'plot_lsline' adds the linear fit line to the data (default = 0);
% set 'color_face' to 0 to remove the shading on the prediction interval
% plots (default = 1);
% Written by Tom Marshall. www.tomrmarshall.com

function [h, ax] = twod_raincloud(X, cl, plot_lsline, color_face)

% defaults
if ~exist('color_face','var')
    color_face = 1;   
end

if ~exist('plot_lsline','var')  
    plot_lsline = 0;  
end

% scatterplot of the two variables
ax(1) = subplot(4, 4, [2 3 4 6 7 8 10 11 12]); hold on
for i = 1:length(X)
    
    h.s{i} = scatter(X{i}(1,:), X{i}(2,:));
    h.s{i}.SizeData         = 32;
    h.s{i}.MarkerFaceColor  = cl(i,:);
    h.s{i}.MarkerFaceAlpha  = 0.5;
    h.s{i}.MarkerEdgeColor  = 'none';
    
end

% least squares lines
if plot_lsline
    
    h.lsl = lsline;
    
    for i = 1:length(X)
        
        h.lsl(end + 1 - i).LineWidth    = 2;
        h.lsl(end + 1 - i).Color        = [cl(i,:) 0.2];
        h.lsl(end + 1 - i).LineStyle    = '--';
        
    end
    
end

% get limits so we can set density plot to same limits
xl = get(gca, 'XLim');
yl = get(gca, 'YLim');

% plot the prediction intervals
for i=1:length(X)
    
    xvals   = linspace(xl(1), xl(2), 10);
    ft      = fit(X{i}(1, :)', X{i}(2, :)', 'poly1');
    predi   = predint(ft, xvals, 0.95, 'functional', 'off');
    h.l{i}  = plot(ax(1), xvals, predi, 'col', [cl(i, :) 0.2], 'LineStyle', '--','LineWidth',2);
    x2      = [xvals fliplr(xvals)];
    predi2  = [predi(:,1); flipud(predi(:,2))];
    
    if color_face
        
        h.f{i} = fill(x2,predi2,cl(i,:));
        h.f{i}.FaceAlpha = 0.05;
        h.f{i}.LineStyle = 'none';
        
    end
    
end

set(gca, 'XLim', xl);
set(gca, 'YLim', yl);
hold off

% x-axis density plots
ax(2) = subplot(4, 4, [14 15 16]);

% calculate kernel densities
for i = 1:length(X)
    
    [f{i}, xi{i}] = ksdensity(X{i}(1,:));
    
    % density plots
    h.dx{i} = area(xi{i}, f{i}); hold on
    set(h.dx{i}, 'FaceColor', cl(i,:));
    set(h.dx{i}, 'EdgeColor', 'none');
    set(h.dx{i}, 'FaceAlpha', 0.5);
    
end

% flip upside down and set limits correctly
set(gca, 'Ydir', 'reverse');
set(gca, 'XLim', xl);
set(gca, 'Visible', 'off');

hold off

% y-axis density plots
ax(3) = subplot(4, 4, [1 5 9]);

% calculate kernel densities
for i = 1:length(X)
    
    [f{i}, xi{i}] = ksdensity(X{i}(2,:));
    
    % density plots
    h.dy{i} = area(xi{i}, f{i}); hold on
    set(h.dy{i}, 'FaceColor', cl(i,:));
    set(h.dy{i}, 'EdgeColor', 'none');
    set(h.dy{i}, 'FaceAlpha', 0.5);
    
end

% rotate this and set the limits correctly
set(gca, 'XLim', yl);
view([-90 -90])
set(gca, 'Xdir', 'reverse');
set(gca, 'Visible', 'off');

hold off

end