function prt_plot_prediction_reg_bar(PRT, model, axes_handle)
% FORMAT prt_plot_prediction_reg_bar(PRT, model, axes_handle)
%
% This function plots the bar plot that appears on prt_ui_results
% Inputs:
%       PRT             - data/design/model structure (it needs to contain
%                         at least one estimated model).
%       model           - the number of the model that will be ploted
%       axes_handle     - (Optional) axes where the plot will be displayed
%
% Output:
%       None
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id: prt_plot_prediction_reg_bar.m 706 2013-06-07 14:33:34Z cphillip $

nfold = length(PRT.model(model).output.fold);

%If no axes_handle is given, create a new window
if ~exist('axes_handle', 'var')
    figure;
    axes_handle = axes;
else
    set(axes_handle, 'XScale','linear');
end

cla(axes_handle, 'reset');
preds1 = [];
preds2 = [];
for f = 1:nfold
    preds1 = [preds1; PRT.model(model).output.fold(f).targets];
    preds2 = [preds2; PRT.model(model).output.fold(f).predictions];
end
bar(axes_handle,[preds1 preds2]);
xlabel(axes_handle,'subjects','FontWeight','bold');
ylabel(axes_handle,'targets and predictions','FontWeight','bold');
legend(axes_handle,{'Target', 'Predicted'});
