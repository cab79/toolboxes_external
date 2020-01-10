function prt_plot_nested_cv(PRT, model, fold, axes_handle)
% FORMAT prt_plot_nested_cv(PRT, model, fold, axes_handle)
%
% Plots the results of the nested cv that appear on prt_ui_results.
%
%
% Inputs:
%       PRT             - data/design/model structure (it needs to contain
%                         at least one estimated model).
%       model           - the number of the model that will be ploted
%       fold            - the number of the fold
%       axes_handle     - (Optional) axes where the plot will be displayed
%
% Output:
%       None
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Matos Monteiro
% $Id$



% Check machine and set the labels an axes
logscale = 0;
switch PRT.model(model).input.machine.function
    case {'prt_machine_svm_bin','prt_machine_sMKL_cla'}
        x_label = 'C';
        y_label = 'Balanced Accuracy (%)';
        
        %If no axes_handle is given, create a new window
        if ~exist('axes_handle', 'var')
            figure;
            axes_handle = axes('XMinorTick','on');
            logscale = 1;
        else
            % Clear EVERYTHING in the UI before defining the axes
            cla(axes_handle, 'reset');
            set(axes_handle,'XMinorTick','on');
            logscale = 1;
        end
        box(axes_handle,'on');
        hold(axes_handle,'all');
        
    case 'prt_machine_sMKL_reg'
        x_label = 'Args';
        y_label = 'MSE';
        
        %If no axes_handle is given, create a new window
        if ~exist('axes_handle', 'var')
            figure;
            axes_handle = axes('XMinorTick','on');
            logscale = 1;
        else
            % Clear EVERYTHING in the UI before defining the axes
            cla(axes_handle, 'reset');
            set(axes_handle, 'XMinorTick','on');
            logscale = 1;
        end
        box(axes_handle,'on');
        hold(axes_handle,'all');
        
        
    case 'prt_machine_krr'
        x_label = 'Args';
        y_label = 'MSE';
        
        %If no axes_handle is given, create a new window
        if ~exist('axes_handle', 'var')
            figure;
            axes_handle = axes;
            logscale = 1;
        else
            % Clear EVERYTHING in the UI before defining the axes
            cla(axes_handle, 'reset');
            logscale = 1;
        end
        
        
        
    case 'prt_machine_wip'
        x_label = 'mu';
        y_label = 'C';
        z_label = 'Balanced Accuracy (%)';
        
        % If no axes_handle is given, create a new window
        if ~exist('axes_handle', 'var')
            figure;
            axes_handle = axes;
        else
            % Clear EVERYTHING in the UI before defining the axes
            cla(axes_handle, 'reset');
            set(axes_handle, 'XScale','linear', 'XMinorTick','on', 'YMinorTick','on');
            logscale = 1;
        end
        
    otherwise
        error('Machine not currently supported for nested CV');
end


cla(axes_handle)
rotate3d off
set(axes_handle,'Color',[1,1,1])
pos=get(axes_handle,'Position');
set(axes_handle,'Position',[pos(1) pos(2) 0.9*pos(3) pos(4)])



% Check if it's a 2 parameter optimisation problem
if strcmp(PRT.model(model).input.machine.function, 'prt_machine_wip')
    
    if fold == 1
        
        nfold = length(PRT.model(model).output.fold);
        
        % Get all function values
        c = unique(PRT.model(model).output.fold(fold).param_effect.param(1,:));
        mu = unique(PRT.model(model).output.fold(fold).param_effect.param(2,:));
        
        for i = 1:nfold
            f(:,:,i) = PRT.model(model).output.fold(i).param_effect.vary_param;
        end
        
        f_mean = mean(f, 3);
        %         f_std = std(f, 0, 3);
        
        f_mean = 100.*f_mean;
        
        % Plot points
        
        %         subplot(2,1,1);
        % TODO: Put Logscale on the Y axis
        axes_handle = image(f_mean, 'CDataMapping', 'scaled', 'XData', mu, 'YData', log10(c));
        % set(axes_handle,'Yscale','log','Ydir','normal');
        axes_color = colorbar;
        title('Mean')
        %         subplot(2,1,2);
        %         axes_handle = image(f_std, 'CDataMapping', 'scaled', 'XData', [min(mu), max(mu)], 'YData', [min(c) max(c)]);
        %         title('Standard Deviation')
        %         colorbar;
        
        
        % Properties
        xlabel(x_label,'FontWeight','bold');
        ylabel(y_label,'FontWeight','bold');
        ylabel(axes_color, z_label,'FontWeight','bold');
        
        
        
        
        % TODO: Try to do it this way instead
        % Include the str information: http://code.izzid.com/2007/08/19/How-to-make-a-3D-plot-with-errorbars-in-matlab.html
        %==================================================================
        %         % TODO: Delete these variables
        %         d_mean = f_mean;
        %         d_std = f_std;
        %
        %         % convert matrices to vectors
        %         f_mean = reshape(f_mean', 1, size(f_mean,1)*size(f_mean,2));
        %         f_std = reshape(f_std', 1, size(f_std,1)*size(f_std,2));
        %
        %         % make mu and x vectors of the same size as f
        %         l_mu = length(mu);
        %         l_c = length(c);
        %         mu = repmat(mu, 1, l_c);
        %         c = repmat(c, l_mu, 1);
        %         c = reshape(c, 1, length(f_mean));
        %
        %         rotate3d on
        %         hold off
        %
        %         axes_handle = plot3(mu, c, f_mean, '.k', 'MarkerSize', 25);
        %         set(axes_handle, 'YScale','log','YMinorTick','on');
        %         %         axes_handle = axes('YScale','log','YMinorTick','on');
        %
        %         hold on
        %         % Draw errorbar for each point
        %         for i = length(f_mean)
        %             c_error = [c(i); c(i)];
        %             mu_error = [mu(i); mu(i)];
        %
        %             f_mean_min = f_mean(i) + f_std(i);
        %             f_mean_max = f_mean(i) - f_std(i);
        %             f_mean_error = [f_mean_min; f_mean_max];
        %
        %             % draw vertical error bar
        %             axes_handle = plot3(mu_error, c_error, f_mean_error, '-k','LineWidth', 2);
        %
        %         end
        %
        %         %         TODO: Finish this!
        %
        
        %==================================================================
        
    else
        
        % Get function values
        c = unique(PRT.model(model).output.fold(fold-1).param_effect.param(1,:));
        mu = unique(PRT.model(model).output.fold(fold-1).param_effect.param(2,:));
        
        f = PRT.model(model).output.fold(fold-1).param_effect.vary_param;
        
        f = 100.*f;
        
        % Plot points
        axes_handle = image(f, 'CDataMapping', 'scaled', 'XData', mu, 'YData', log10(c));
        axes_color = colorbar;
        
        % Properties
        xlabel(x_label,'FontWeight','bold');
        ylabel(y_label,'FontWeight','bold');
        ylabel(axes_color, z_label,'FontWeight','bold');
        
        
    end
    
    
else % It's a 1 parameter optimisation problem
    
    
    if fold == 1
        
        nfold = length(PRT.model(model).output.fold);
        
        % Get function values
        x = PRT.model(model).output.fold(fold).param_effect.param;
        f = zeros(nfold, length(x));
        
        % Get mean f values
        for i = 1:nfold
            f(i,:) = PRT.model(model).output.fold(i).param_effect.vary_param;
            % Get the chosen optimal values
            x_opt(i) = PRT.model(model).output.fold(i).param_effect.opt_param;
        end
        
        if strcmp(PRT.model(model).input.type, 'classification')
            f = 100.*f; % Convert to percentage
        end
        f_mean = mean(f);
        f_std = std(f);
        
        % get frequencies of optimal values
        x_opt = hist(x_opt, x)./size(f,1);
        
        
         % Plot
        if logscale
            x = log10(x);
        end
        
        % general properties of the plots
        markersize = 10;
        f_min = 0;
        f_max = 108;

        hold on
        [hax,hbar,hline] = plotyy(x,x_opt*100,x,mean(f),'bar','plot');
        errorbar(axes_handle, x, f_mean, f_std, '.k', 'linewidth', 2);
        set(hbar,'BarWidth',0.5,'FaceColor',[0.5 0.8 0.5])
        set(hline,'Color','k','Linewidth',1)
        set(hax(1),'YColor',[0.1,0.6,0.1])
        set(hax(2),'YColor',[0,0,0])
        hold off
        
        % Properties
       
        ylabel(hax(2), y_label,'FontWeight','bold');
        ylabel(hax(1),'Frequency of selection (%)','FontWeight','bold');
        if logscale 
             xlabel(axes_handle, [x_label, ' (log 10)'],'FontWeight','bold');
        else
             xlabel(axes_handle, x_label,'FontWeight','bold');
        end
        axis(hax(1), [min(x)-0.2*abs(min(x)) max(x)+0.2*abs(max(x)) f_min f_max]);
        axis(hax(2), [min(x)-0.2*abs(min(x)) max(x)+0.2*abs(max(x)) f_min f_max]);
        set(hax(2),'XTickLabel',{})      
        set(hax(2),'YTickLabel',{})
        a=get(hax(1),'YTick');
        b=get(hax(1),'YTickLabel');
        set(hax(2),'YTick',a);
        set(hax(2),'YTickLabel',b);         
        
    else
        
        % Get all function values
        x = PRT.model(model).output.fold(fold-1).param_effect.param;
        f = PRT.model(model).output.fold(fold-1).param_effect.vary_param;
        if strcmp(PRT.model(model).input.type, 'classification')
            f = 100.*f; % Convert to percentage
        end
        
        % Get optimal function values
        switch PRT.model(model).input.type
            case 'classification'
                x_opt = find(f==max(f));
            case 'regression'
                x_opt = find(f==min(f));
            otherwise
                error('Type of model not recognised');
        end
        
        % general properties of the plots
        markersize = 10;
        switch PRT.model(model).input.type
            case 'classification'
                f_min = 0;
                f_max = 108;
            case 'regression'
                f_min = min(f(:))-0.1*min(f(:));
                f_max = max(f(:))+0.1*min(f(:));
            otherwise
                error('Type of model not recognised');
        end
        
        if logscale
            x = log10(x);
        end
        
        % Plot all points
        hold on
        plot(axes_handle, x, f, '-xk', 'markersize', markersize, 'linewidth', 1);
        % Plot the optimal on top of the original
        opt_handle = plot(axes_handle, x(x_opt), f(x_opt), 'xr', 'markersize', markersize, 'linewidth', 3);
        hold off
        
        % Properties
        if logscale
            x_label = [x_label,' (log 10)'];
        end
        xlabel(axes_handle, x_label,'FontWeight','bold');
        ylabel(axes_handle, y_label,'FontWeight','bold');
        legend(opt_handle, 'Optimal value(s)');
        axis(axes_handle, [min(x) max(x) f_min f_max]);
        
    end
    
end

end
