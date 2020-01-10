function prt_plot_histograms(PRT, model, fold, axes_handle)
% FORMAT prt_plot_histograms(PRT, model, fold, axes_handle)
%
% This function plots the histogram that appears on prt_ui_results.
%
% The maximum number of classes that can be ploted is 7. However, this can
% be increased by editing the function. Just add more colours to the
% colourList variable.
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

% Written by M. J. Rosa
% $Id: prt_plot_histograms.m 706 2013-06-07 14:33:34Z cphillip $

nfold = length(PRT.model(model).output.fold);


%Check the number of classes
nClasses = length(PRT.model(model).input.class);

if fold == 1
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
    
    
    for i=1:nClasses
        targval(:, i) = targets == i;
    end
    
else
    % if folds wise
    targets = PRT.model(model).output.fold(fold-1).targets;
    
    for i=1:nClasses
        targval(:, i) = targets == i;
    end
    
    
    if isfield(PRT.model(model).output.fold(fold-1),'func_val')
        fVals  = PRT.model(model).output.fold(fold-1).func_val;
        fVvals_exist = 1;
    else
        fVvals_exist = 0;
        fVals  = PRT.model(model).output.fold(fold-1).predictions;
    end
end


%Make list of classes
for i=1:nClasses
    classNames{i} = PRT.model(model).input.class(i).class_name;
end
%If you want to use more classes, just add more colours to the list bellow
colourList = {'black','red', 'blue', 'green', 'cyan', 'magenta', 'yellow'};


%If no axes_handle is given, create a new window
if ~exist('axes_handle', 'var')
    figure;
    axes_handle = axes;
else
    set(axes_handle, 'XScale','linear');
end

cla(axes_handle, 'reset');
rotate3d off
%                 axis xy
set(axes_handle,'Color',[1,1,1])

if nClasses <= length(colourList)
    if fVvals_exist
        classes_used = [];
        for cl=1:nClasses
            func_vals = fVals(targval(:,cl));
            if ~isempty(func_vals)
                if exist('ksdensity','file')==2
                    [f,x] = ksdensity(func_vals,'width',[]);
                    plot(axes_handle,x,f,colourList{cl},'LineWidth',2);
                    hold(axes_handle,'on')
                else
                    % can't plot density, be happy with a histogram
                    [myHist,myX]=hist(func_vals,100);
                    bar(axes_handle,myX,myHist,colourList{cl});
                    hold(axes_handle,'on')
                end
                if cl == nClasses, hold(axes_handle,'off'); end
                classes_used = [classes_used,cl]; % makes a list of the classes used in the plot
            end
        end
        
        %Make list of classes used to put in the legend
        nClasses_used = length(classes_used);
        for i=1:nClasses_used
            classNames_legend{i} = classNames{classes_used(i)};
        end
        
        legend(axes_handle,classNames_legend);
        
        xlabel(axes_handle,'function value','FontWeight','bold');
    else
        % do nothing, no func_val available
    end
    
else
    error(['Too many classes, Max number of classes: ' num2str(length(colourList))]);
end
