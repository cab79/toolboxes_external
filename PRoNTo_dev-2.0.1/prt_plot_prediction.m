function prt_plot_prediction(PRT, model, fold, marker_size, axes_handle)
% FORMAT prt_plot_prediction(PRT, model, fold, marker_size, axes_handle)
%
% This function plots the prediction plot that appears on prt_ui_results
% Inputs:
%       PRT             - data/design/model structure (it needs to contain
%                         at least one estimated model).
%       model           - the number of the model that will be ploted
%       fold            - the number of the fold
%       marker_size     - (Optional) the size of the markers in the plot,
%                         the default is 7
%       axes_handle     - (Optional) axes where the plot will be displayed
%
% Output:
%       None
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id: prt_plot_prediction.m 706 2013-06-07 14:33:34Z cphillip $


nfold = length(PRT.model(model).output.fold);
bcl = 0;


fVals   = [];
targets = [];

for f = 1:nfold,
    targets = [targets;PRT.model(model).output.fold(f).targets];
    if isfield(PRT.model(model).output.fold(f),'func_val')
        fVvals_exist = 1;
        fVals  = [fVals;PRT.model(model).output.fold(f).func_val];
    else
        fVvals_exist = 0;
        fVals  = [fVals;...
            PRT.model(model).output.fold(f).predictions];
    end
end
tarmax = max(targets);
if fold==1
    targpos = targets == 1;
    if length(find(targpos))~=length(targets) && ~isempty(find(targpos))
        %both classes are present across all folds (for class labels)
        bcl = 1;
    end
end
    
if fold>1
    % if folds wise
    targets = PRT.model(model).output.fold(fold-1).targets;
    targpos = targets == tarmax;
    if isfield(PRT.model(model).output.fold(fold-1),'func_val')
        fVals  = PRT.model(model).output.fold(fold-1).func_val;
        fVvals_exist = 1;
    else
        fVvals_exist = 0;
        fVals  = PRT.model(model).output.fold(fold-1).predictions;
    end
end


%Defined the marker size, if no value is given
if ~exist('marker_size', 'var')
    marker_size = 7;
end
%If no axes_handle is given, create a new window
if ~exist('axes_handle', 'var')
    figure;
    axes_handle = axes;
else
    set(axes_handle, 'XScale','linear');
end

cla(axes_handle, 'reset');
rotate3d off
colorbar('peer',axes_handle,'off')
set(axes_handle,'Color',[1,1,1])
% predictions
if fVvals_exist
    if fold == 1
        foldlabels = 1:nfold;
        for f = 2:nfold+1
            targets = PRT.model(model).output.fold(f-1).targets;
            targpos = targets == tarmax; %both for SVM and GP, linear
            fVals   = PRT.model(model).output.fold(f-1).func_val;
            func_valsc1 = fVals(targpos);
            func_valsc2 = fVals(~targpos);
            yc1 = (f-1)*ones(length(func_valsc1),1);
            yc2 = (f-1)*ones(length(func_valsc2),1);
            if f==2
                maxfv = max(abs([func_valsc1;func_valsc2]));
            else
                maxtmp = max(abs([func_valsc1;func_valsc2]));
                if maxfv < maxtmp, maxfv = maxtmp; end
            end
            pl1 = plot(axes_handle,func_valsc1,yc1,'kx','MarkerSize',marker_size);
            hold(axes_handle,'on');
            if ~isempty(yc1), isyc1 = 1; plot1 = pl1; else isyc1 = 0; end
            pl2 = plot(axes_handle,func_valsc2,yc2,'ro','MarkerSize',marker_size);
            hold(axes_handle,'on');
            if ~isempty(yc2), isyc2 = 1; plot2 = pl2; else isyc2 = 0; end
        end
    else
        foldlabels  = fold-1;
        func_valsc1 = fVals(targpos);
        func_valsc2 = fVals(~targpos);
        yc1 = (fold-1)*ones(length(func_valsc1),1);
        yc2 = (fold-1)*ones(length(func_valsc2),1);
        maxfv = max(abs([func_valsc1;func_valsc2]));
        pl1 = plot(axes_handle,func_valsc1,yc1,'kx','MarkerSize',marker_size);
        hold(axes_handle,'on');
        if ~isempty(yc1), isyc1 = 1; plot1 = pl1; else isyc1 = 0; end
        pl2 = plot(axes_handle,func_valsc2,yc2,'ro','MarkerSize',marker_size);
        hold(axes_handle,'on');
        if ~isempty(yc2), isyc2 = 1; plot2 = pl2; else isyc2 = 0;  end
    end
    % Change the x axis for gaussian process or RT - change in
    % the future
    y = [0:nfold+1]';
    if strcmp(PRT.model(model).input.machine.function,'prt_machine_gpml');
        x = 0.5*ones(nfold+2,1);
        plot(axes_handle,x,y,'--','Color',[1 1 1]*.6);
        xlim(axes_handle,[0 1]);
    elseif strcmp(PRT.model(model).input.machine.function,'prt_machine_RT_bin')
        % nothing to do - just leave auto scaling
    else
        x = zeros(nfold+2,1);
        plot(axes_handle,x,y,'--','Color',[1 1 1]*.6);
        xlim(axes_handle,[-maxfv-0.5 maxfv+0.5]);
    end
    ylim(axes_handle,[0 nfold+1.3]);
    xlabel(axes_handle,'function value','FontWeight','bold');
    h=ylabel(axes_handle,'fold','FontWeight','bold');
    set(h,'Rotation',90)
    
    %These 2 lines were added. See if it's better to give the class names
    %as input
    classNames{1} = PRT.model(model).input.class(2).class_name;
    classNames{2} = PRT.model(model).input.class(1).class_name;
    
    
    if (isyc1 && isyc2) || bcl
        legend([plot1,plot2],classNames,'Color',[1,1,1]);
    else
        if isyc1
            legend(plot1,classNames{1});
        else
            legend(plot2,classNames{2});
        end
    end
    set(axes_handle,'YTick',foldlabels)
    hold(axes_handle,'off');
    set(axes_handle,'Color',[1,1,1],'Visible','on')
    %                 axis normal
    %                 axis xy
    title(axes_handle,'')
else
    set(axes_handle,'Color',[1,1,1])
    beep
    disp('No function values to display!')
end
