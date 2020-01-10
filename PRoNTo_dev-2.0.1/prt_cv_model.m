function [outfile,PRT]=prt_cv_model(PRT,in)
% Function to run a cross-validation structure on a given model
%
% Inputs:
% -------
% PRT:             data structure
% in.fname:        filename for PRT.mat (string)
% in.model_name:   name for this model (string)
%
% Outputs:
% --------
% Writes the following fields in the PRT data structure:
%
% PRT.model(m).output.fold(i).targets:     targets for fold(i)
% PRT.model(m).output.fold(i).predictions: predictions for fold(i)
% PRT.model(m).output.fold(i).stats:       statistics for fold(i)
% PRT.model(m).output.fold(i).{custom}:    optional fields
%
% Notes:
% ------
% The PRT.model(m).input fields are set by prt_init_model, not by
% this function
%
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

prt_dir = char(regexprep(in.fname,'PRT.mat', ''));

% Get index of specified model
mid = prt_init_model(PRT, in);

% configure some variables
CV       = PRT.model(mid).input.cv_mat;     % CV matrix
n_folds  = size(CV,2);                      % number of CV folds

% targets
if isfield(PRT.model(mid).input,'include_allscans') && ...
        PRT.model(mid).input.include_allscans
    t = PRT.model(mid).input.targ_allscans;
    % Get covariates if GLM required
    if any(ismember(PRT.model(mid).input.operations,5))
        cov = PRT.model(mid).input.cov_allscans;
    else
        cov=[];
    end
else
    t = PRT.model(mid).input.targets;
    % Get covariates if GLM required
    if any(ismember(PRT.model(mid).input.operations,5))
        cov = PRT.model(mid).input.covar;
    else
        cov=[];
    end
end

%get number of classes
if strcmpi(PRT.model(mid).input.type,'classification')
    nc=max(unique(t));
else
    nc=[];
end
fdata.nc = nc;

%load kernels and get the used sample in this model
[Phi_all,ID] = prt_getKernelModel(PRT,prt_dir,mid);


% Begin cross-validation loop
% -------------------------------------------------------------------------
PRT.model(mid).output=struct();
PRT.model(mid).output.fold = struct();
for f = 1:n_folds
    disp ([' > running CV fold: ',num2str(f),' of ',num2str(n_folds),' ...'])
    % configure data structure for prt_cv_fold
    fdata.ID      = ID;
    fdata.mid     = mid; %index of model
    fdata.CV      = CV(:,f);
    fdata.Phi_all = Phi_all; %kernel(s)
    fdata.t       = t; %targets
    if ~isempty(cov)
        fdata.cov = cov;
    end
    
    % Nested CV for hyper-parameter optimisation or feature selection
    if isfield(PRT.model(mid).input,'use_nested_cv')
        if PRT.model(mid).input.use_nested_cv
            if f==1 && isempty(PRT.model(mid).input.nested_param)
                beep
                warning('No parameter range specified for optimization, using defaults.')
            end
            [out] = prt_nested_cv(PRT, fdata);
            PRT.model(mid).output.fold(f).param_effect = out;
            PRT.model(mid).input.machine.args = out.opt_param;
        end
    end
    
    % compute the model for this CV fold
    [model, targets] = prt_cv_fold(PRT,fdata);
    
    %for classification check that for each fold, the test targets have been trained
    if strcmpi(PRT.model(mid).input.type,'classification')
        if ~all(ismember(unique(targets.test),unique(targets.train)))
            beep
            disp('At least one class is in the test set but not in the training set')
            disp('Abandoning modelling, please correct class selection/cross-validation')
            return
        end
    end
    
    % compute stats
    stats = prt_stats(model, targets.test, nc); %targets.train
    
    % update PRT
    PRT.model(mid).output.fold(f).targets     = targets.test;
    PRT.model(mid).output.fold(f).predictions = model.predictions(:);
    PRT.model(mid).output.fold(f).stats       = stats;
    % copy other fields from the model
    flds = fieldnames(model);
    for fld = 1:length(flds)
        fldnm = char(flds(fld));
        if ~strcmpi(fldnm,'predictions')
            PRT.model(mid).output.fold(f).(fldnm)=model.(fldnm);
        end
    end
end


% Model level statistics (across folds)
ttt             = vertcat(PRT.model(mid).output.fold(:).targets);
m.type        = PRT.model(mid).output.fold(1).type;
m.predictions = vertcat(PRT.model(mid).output.fold(:).predictions);
%m.func_val    = [PRT.model(mid).output.fold(:).func_val];
stats         = prt_stats(m,ttt(:),nc);

PRT.model(mid).output.stats=stats;


% Save PRT containing machine output
% -------------------------------------------------------------------------
outfile = [prt_dir, filesep,'PRT.mat'];
disp('Updating PRT.mat.......>>')
if spm_check_version('MATLAB','7') < 0
    save(outfile,'-V6','PRT');
else
    save(outfile,'PRT');
end
end


