function [PRT, CV, ID] = prt_model(PRT,in) % CAB line 154
% Function to configure and build the PRT.model data structure
%
% Input:
% ------
%   PRT fields:
%   model.fs(f).fs_name:     feature set(s) this CV approach is defined for
%   model.fs(f).fs_features: feature selection mode ('all' or 'mask')
%   model.fs(f).mask_file:   mask for this feature set (fs_features='mask')
%
%   in.fname:      filename for PRT.mat
%   in.model_name: name for this cross-validation structure
%   in.type:       'classification' or 'regression'
%   in.use_kernel: does this model use kernels or features?
%   in.operations: operations to apply before prediction
%
%   in.fs(f).fs_name:     feature set(s) this CV approach is defined for
%
%   in.class(c).class_name
%   in.class(c).group(g).subj(s).num
%   in.class(c).group(g).subj(s).modality(m).mod_name
%   EITHER: in.class(c).group(g).subj(s).modality(m).conds(c).cond_name
%   OR:     in.class(c).group(g).subj(s).modality(m).all_scans
%   OR:     in.class(c).group(g).subj(s).modality(m).all_cond
%
%   in.cv.type:     type of cross-validation ('loso','losgo','custom')
%   in.cv.mat_file: file specifying CV matrix (if type='custom');
%
% Output:
% -------
%
%   This function performs the following functions:
%      1. populates basic fields in PRT.model(m).input
%      2. computes PRT.model(m).input.targets based on in.class(c)...
%      3. computes PRT.model(m).input.samp_idx based on targets
%      4. computes PRT.model(m).input.cv_mat based on the labels and CV spec
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

% Populate basic fields in PRT.mat
% -------------------------------------------------------------------------
[modelid, PRT] = prt_init_model(PRT,in);

% specify model type and feature sets
PRT.model(modelid).input.type = in.type;
if strcmp(in.type,'classification')
    for c = 1:length(in.class)
        PRT.model(modelid).input.class(c) = in.class(c);
    end
end

for f = 1:length(in.fs)
    fid = prt_init_fs(PRT,in.fs(f));
    
    if length(PRT.fs(fid).modality) > 1 && length(in.fs) > 1
        error('prt_model:multipleFeatureSetsAppliedAsSamplesAndAsFeatures',...
            ['Feature set ',in.fs(f).fs_name,' contains multiple modalities ',...
            'and job specifies that multiple feature sets should be ',...
            'supplied to the machine. This usage is not supported.']);
    end
    
    PRT.model(modelid).input.fs(f).fs_name = in.fs(f).fs_name;
end

% compute targets and samp_idx
% -------------------------------------------------------------------------
if strcmp(in.type,'classification')
    [targets, samp_idx, t_allscans, samp_allscans,covar,cov_all] = compute_targets(PRT, in);
else
    [targets, samp_idx, t_allscans,covar,cov_all] = compute_target_reg(PRT, in);
end
%[afm]
if isfield(in,'include_allscans') && in.include_allscans   
    PRT.model(modelid).input.samp_idx = samp_allscans;
    PRT.model(modelid).input.include_allscans = in.include_allscans;
else
    PRT.model(modelid).input.samp_idx = samp_idx;
    PRT.model(modelid).input.include_allscans = false;
end
PRT.model(modelid).input.targets          = targets;
PRT.model(modelid).input.targ_allscans    = t_allscans;
PRT.model(modelid).input.covar            = covar;
PRT.model(modelid).input.cov_allscans     = cov_all;

% compute cross-validation matrix and specify operations to apply
% -------------------------------------------------------------------------
if isfield(in.cv,'k')
    PRT.model(modelid).input.cv_k=in.cv.k;
else
    PRT.model(modelid).input.cv_k = 0;
end  
[CV,ID] = prt_compute_cv_mat(PRT,in, modelid);
PRT.model(modelid).input.cv_mat     = CV;
PRT.model(modelid).input.cv_type=in.cv.type;
% Deal with nested CV parameters
if isfield(in.cv,'type_nested') && ~isempty(in.cv.type_nested)
    PRT.model(modelid).input.cv_type_nested = in.cv.type_nested;
end
if isfield(in.cv,'k_nested') && ~isempty(in.cv.k_nested)
    PRT.model(modelid).input.cv_k_nested = in.cv.k_nested;
end
if isfield(in.cv,'nested_param') && ~isempty(in.cv.nested_param)
    PRT.model(modelid).input.nested_param = in.cv.nested_param;
end

PRT.model(modelid).input.operations = in.operations;

% Save PRT.mat
% -------------------------------------------------------------------------
disp('Updating PRT.mat.......>>')
if spm_check_version('MATLAB','7') >= 0
    save(in.fname,'-V7','PRT');
else
    save(in.fname,'-V6','PRT');
end

end

%% -------------------------------------------------------------------------
% Private Functions
% -------------------------------------------------------------------------

function [targets, samp_idx, t_all samp_all,covar,cov_all] = compute_targets(PRT, in)
% Function to compute the prediction targets. Also does some error checking

% Set the reference feature set
fid = prt_init_fs(PRT, in.fs(1));
ID  = PRT.fs(fid).id_mat;
n   = size(ID,1);

% Check the feature sets have the same number of samples (eg for MKL).
if length(in.fs) > 1
    for f = 1:length(in.fs)
        fid = prt_init_fs(PRT, in.fs(f));
        if size(PRT.fs(fid).id_mat,1) ~= n
            error('prt_model:sizeOfFeatureSetsDiffer',...
                ['Multiple feature sets included, but they have different ',...
                'numbers of samples']);
        end
    end
end

modalities = {PRT.masks(:).mod_name};
groups     = {PRT.group(:).gr_name};

t_all    = zeros(n,1);
samp_all = zeros(n,1);
cov_all = zeros(n,1);

% CAB - assumes binary, not multiclass
class_labels = repmat([1 2],1,length(in.class)/2);

for c = 1:length(in.class)
    
    % CAB
    cl = class_labels(c);
    
    % groups
    for g = 1:length(in.class(c).group)
        gr_name = in.class(c).group(g).gr_name;
        if any(strcmpi(gr_name,groups))
            gid = find(strcmpi(gr_name,groups));
        else
            error('prt_model:groupNotFoundInPRT',...
                ['Group ',gr_name,' not found in PRT.mat']);
        end
        
        % subjects
        for s = 1:length(in.class(c).group(g).subj)
            sid = in.class(c).group(g).subj(s).num;
            % modalities
            for m = 1:length(in.class(c).group(g).subj(s).modality)
                mod_name = in.class(c).group(g).subj(s).modality(m).mod_name;
                if any(strcmpi(mod_name,modalities))
                    mid = find(strcmpi(mod_name,modalities));
                else
                    error('prt_model:groupNotFoundInPRT',...
                        ['Modality ',mod_name,' not found in PRT.mat']);
                end
                
                if isfield(in.class(c).group(g).subj(s).modality(m), 'all_scans')
                    % check whether this was included in the feature set
                    % using 'all conditions' (which is invalid)
                    if strcmpi(PRT.fs(fid).modality(m).mode,'all_cond')
                        error('prt_model:fsIsAllCondModelisAllScans',...
                            ['''All scans'' selected for subject ',num2str(s),...
                            ', group ',num2str(g), ', modality ', num2str(m),...
                            ' but the feature set was constructed using ',...
                            '''All conditions''. This syntax is invalid. ',...
                            'Please use ''All Conditions'' instead.']);
                    end
                    
                    % otherwise add all scans for each subject
                    %[afm] idx = ID(:,1) == gid & ID(:,2) == s & ID(:,3) == mid;
                    idx = ID(:,1) == gid & ID(:,2) == sid & ID(:,3) == mid;
                    t_all(idx) = cl;
                    if any(ismember(in.operations, 5)) %Get covariates
                        cov_all(idx) = PRT.group(gid).subject(sid).modality(mid).covar;
                    end
                else % conditions have been specified
                    % check whether conditions were specified in the design
                    if ~isfield(PRT.group(gid).subject(sid).modality(mid).design,'conds')
                        error('prt_model:conditionsSpecifiedButNoneInDesign',...
                            ['Conditions selected for subject ',num2str(s),...
                            ', class ',num2str(c),', group ',num2str(g), ...
                            ', modality ', num2str(m),' but there are none in the design. ',...
                            'Please use ''All Scans'' or adjust design.']);
                    end
                    %[afm]sid = in.class(c).group(g).subj(s).num;
                    conds     = {PRT.group(gid).subject(sid).modality(mid).design.conds(:).cond_name};
                    
                    
                    if isfield(in.class(c).group(g).subj(s).modality(m), 'all_cond')
                        % all conditions
                        for cid = 1:length(conds)
                            idx = ID(:,1) == gid & ID(:,2) == sid & ID(:,3) == mid & ID(:,4) == cid;
                            t_all(idx) = cl;
                        end
                    else % loop over conditions
                        for cond = 1:length(in.class(c).group(g).subj(s).modality(m).conds)
                            cond_name = in.class(c).group(g).subj(s).modality(m).conds(cond).cond_name;
                            
                            if any(strcmpi(cond_name,conds))
                                cid = find(strcmpi(cond_name,conds));
                            else
                                error('prt_model:groupNotFoundInPRT',...
                                    ['Condition ',cond_name,' not found in PRT.mat']);
                            end
                            
                            idx = ID(:,1) == gid & ID(:,2) == sid & ID(:,3) == mid & ID(:,4) == cid;
                            t_all(idx) = cl;
                        end
                    end
                    s_idx_mod = ID(:,1) == gid & ID(:,2) == sid & ID(:,3) == mid;
                    samp_all(s_idx_mod) = 1;
                end
            end
        end
    end
end

samp_idx = find(t_all);
samp_all = find(samp_all);
targets  = t_all(samp_idx);
covar = cov_all(samp_idx);
end


function [targets, samp_idx, targ_allscans,covar,cov_all]=compute_target_reg(PRT, in)
% Function to compute the prediction targets. Not much error checking yet

% Set the reference feature set
fid = prt_init_fs(PRT, in.fs(1));
ID  = PRT.fs(fid).id_mat;
n   = size(ID,1);

modalities = {PRT.masks(:).mod_name};
groups     = {PRT.group(:).gr_name};
%t_all = zeros(n,1);
targ_allscans=zeros(n,1);
cov_all = zeros(n,1);
samp_idx=[];
targ_g=[];
covar = [];
for g = 1:length(in.group)
    gr_name = in.group(g).gr_name;
    if any(strcmpi(gr_name,groups))
        gid = find(strcmpi(gr_name,groups));
    else
        error('prt_model:groupNotFoundInPRT',...
            ['Group ',gr_name,' not found in PRT.mat']);
    end
%     nmod=length(in.group(g).subj(1).modality);
    targets=zeros(1,length(in.group(g).subj)); %replace by nmod for multiple targets per subject
    cov=zeros(1,length(in.group(g).subj));
    % subjects
    for s = 1:length(in.group(g).subj)
        %modalities
        for m = 1:length(in.group(g).subj(s).modality)
            mod_name = in.group(g).subj(s).modality(m).mod_name;
            if any(strcmpi(mod_name,modalities))
                mid = find(strcmpi(mod_name,modalities));
            else
                error('prt_model:groupNotFoundInPRT',...
                    ['Modality ',mod_name,' not found in PRT.mat']);
            end
            if m==1 %only one regression target per subject, whatever the number of modalities
                idx = in.group(g).subj(s).num;
                if ~isempty(PRT.group(gid).subject(idx).modality(mid).rt_subj)
                    targets(m,s) = PRT.group(gid).subject(idx).modality(mid).rt_subj;
                else
                    error('prt_model:NoRegressionTarget','No regression target found, correct');
                end
                samp_idx=[samp_idx; find(ID(:,1) == gid & ID(:,2) == idx & ID(:,3) == mid)];
                if any(ismember(in.operations, 5)) %Get covariates
                    cov(m,s) = PRT.group(gid).subject(idx).modality(mid).covar;
                end
            end
        end        
    end
    targ_g=[targ_g;targets(:)];
    covar = [covar;cov(:)];
end
targ_allscans(samp_idx)=targ_g;
targets=targ_g;
cov_all(samp_idx)=covar;
end
