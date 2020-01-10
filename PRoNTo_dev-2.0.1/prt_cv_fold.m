function [model, targets] = prt_cv_fold(PRT, in)
% Function to run a single cross-validation fold 
%
% Inputs:
% -------
% PRT:           data structure
% in.mid:        index to the model we are working on
% in.ID:         ID matrix
% in.CV:         Cross-validation matrix (current fold only)
% in.Phi_all:    Cell array of data matri(ces) (training and test)
% in.t           prediction targets
%
% Outputs:
% --------
% model:         the model returned by the machine
% targets.train: training targets
% targets.test:  test targets
%
% Notes: 
% ------
% The training and test targets output byt this function are not
% necessarily equivalent to the targets that are supplied to the function.
% e.g. some data operations can modify the number of samples (e.g. sample
% averaging). In such cases size(targets.train) ~= size(in.t)
%
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand 
% $Id$

tr_idx = in.CV == 1;
te_idx = in.CV == 2;

[Phi_tr, Phi_te, Phi_tt] = ...
    split_data(in.Phi_all, tr_idx, te_idx, PRT.model(in.mid).input.use_kernel);

% Assemble data structure to supply to machine
cvdata.train      = Phi_tr;
cvdata.test       = Phi_te;
if PRT.model(in.mid).input.use_kernel
    cvdata.testcov    = Phi_tt;
end

% configure basic CV parameters
cvdata.tr_targets = in.t(tr_idx,:);
cvdata.te_targets = in.t(te_idx,:);
cvdata.tr_id      = in.ID(tr_idx,:);
cvdata.te_id      = in.ID(te_idx,:);
cvdata.use_kernel = PRT.model(in.mid).input.use_kernel;
cvdata.pred_type  = PRT.model(in.mid).input.type;

% configure additional CV parameters (e.g. needed to compute a GLM)
cvdata.tr_param = prt_cv_opt_param(PRT, in.ID(tr_idx,:), in.mid);
cvdata.te_param = prt_cv_opt_param(PRT, in.ID(te_idx,:), in.mid);

% Apply any operations specified
ops = PRT.model(in.mid).input.operations(PRT.model(in.mid).input.operations ~=0 );
if any(ismember(ops,5))
    cvdata.tr_cov = in.cov(tr_idx,:);
    cvdata.te_cov = in.cov(te_idx,:);
    posglm = find(ops==5);
    if posglm~=1 % GLM should be first
        idxops = 1:length(ops);
        newidx = setdiff(idxops,posglm);
        ops=[5,ops(newidx)];
    end
end
for o = 1:length(ops)    
    cvdata = prt_apply_operation(PRT, cvdata, ops(o));
end

% train the prediction model
try
    model = prt_machine(cvdata, PRT.model(in.mid).input.machine);
catch err
    warning('prt_cv_fold:modelDidNotReturn',...
        'Prediction method did not return [%s]',err.message);
    model.predictions = zeros(size(cvdata.te_targets));
end

% check that it produced a predictions field
if ~any(strcmpi(fieldnames(model),'predictions'))
    error(['prt_cv_model:machineDoesNotGivePredictions',...
        'Machine did not produce a predictions field']);
end

% does the model alter the target vector (e.g. change its dimension) ?
if isfield(model,'te_targets')
    targets.test = model.te_targets(:);
else
    targets.test = cvdata.te_targets(:);
end
if isfield(model,'tr_targets')
    targets.train = model.tr_targets(:);
else
    targets.train= cvdata.tr_targets(:);
end

end

% -------------------------------------------------------------------------
% Private functions
% -------------------------------------------------------------------------
        
function [Phi_tr Phi_te Phi_tt] = split_data(Phi_all, tr_idx, te_idx, usebf)
% function to split the data matrix into training and test

n_mat = length(Phi_all);

% training
Phi_tr = cell(1,n_mat);
for i = 1:n_mat;
    if usebf
        cols_tr = tr_idx;
    else
        cols_tr = size(Phi_all{i},2);
    end
    
    Phi_tr{i} = Phi_all{i}(tr_idx,cols_tr);

end

% test
Phi_te  = cell(1,n_mat);
Phi_tt = cell(1,n_mat);
if usebf
    cols_tr = tr_idx;
    cols_te = te_idx;
else
    cols_tr = size(Phi_all{i},2);
    %cols_te = size(Phi_all{i},2);
end

for i = 1:length(Phi_all)
    Phi_te{i} = Phi_all{i}(te_idx, cols_tr);
    if usebf
        Phi_tt{i} = Phi_all{i}(te_idx, cols_te);
    else
        Phi_tt{i} = [];
    end
end
end
