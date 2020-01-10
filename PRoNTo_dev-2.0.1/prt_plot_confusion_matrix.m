function prt_plot_confusion_matrix(PRT, model, fold, axes_handle)
% FORMAT prt_plot_confusion_matrix(PRT, model, fold, axes_handle)
%
% This function plots the confusion matrix that appears on prt_ui_results
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
% $Id: prt_plot_confusion_matrix.m 706 2013-06-07 14:33:34Z cphillip $

%If no axes_handle is given, create a new window
if ~exist('axes_handle', 'var')
    figure;
    axes_handle = axes;
else
    set(axes_handle, 'XScale','linear');
end


% confusion matrix
cla(axes_handle, 'reset');
if fold == 1
    mconmat(:,:) = PRT.model(model).output.stats.con_mat;
else
    mconmat(:,:) = PRT.model(model).output.fold(fold-1).stats.con_mat;
end
myH=bar3(axes_handle,mconmat,'detached','w');
nclass = size(mconmat,1);
for j = 1:nclass,
    conLabels{j} = num2str(j);
end
rotate3d on
if fold == 1
    title(axes_handle,sprintf('Confusion matrix: all folds'),'FontWeight','bold');
else
    title(axes_handle,sprintf('Confusion matrix: fold %d',fold-1),'FontWeight','bold');
end
xlabel(axes_handle,'True','FontWeight','bold');
ylabel(axes_handle,'Predicted','FontWeight','bold');
set(axes_handle,'XTick',1:nclass);
set(axes_handle,'XTickLabel',conLabels);
set(axes_handle,'YTick',1:nclass);
set(axes_handle,'YTickLabel',conLabels);
grid(axes_handle,'on');
set(axes_handle,'Color',[0.8 0.8 0.8]);
axis square; axis vis3d; axis tight;
% add values
for foo_row=1:size(mconmat,1)
    for foo_col=1:size(mconmat,2)
        foo_zval=mconmat(foo_row,foo_col);
        if foo_row==foo_col, foo_color='g'; else foo_color='r';end
        text(foo_col,foo_row,foo_zval,num2str(foo_zval),...
            'Color',foo_color);
    end
end
