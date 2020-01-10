function [outfile,PRT] = prt_permutation(PRT, n_perm, modelid, path, flag)
% Function to compute permutation test
%
% Inputs:
% -------
% PRT:     PRT structured including model
% n_perm:  number of permutations
% modelid: model ID
% path:    path
% flag:    boolean variable. set to 1 to save the weights for each
%          permutation. default: 0
%
% Outputs:
% --------
%
% for classification
% permutation.c_acc:        Permuted accuracy per class
% permutation.b_acc:        Permuted balanced accuracy
% permutation.pvalue_b_acc: p-value for c_acc
% permutation.pvalue_c_acc: p-value for b_acc
%
% for regression
% permutation.corr: Permuted correlation
% permutation.mse:  Permuted mean square error
% permutation.pval_corr: p-value for corr
% permutation.pval_r2: p-value for r2;
% permutation.pval_mse:  p-value for mse
% permutation.pval_nmse:  p-value for nmse
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Mourao-Miranda
% $Id$

prt_dir = path;
def_par = prt_get_defaults('paral');
if nargin<5
    flag=0;
end

% % prt_dir = char(regexprep(in.fname,'PRT.mat', ''));

if ~isfield(PRT,'model')
    beep
    disp('No model found in this PRT.mat');
    return
else
    if ~isfield(PRT.model,'output')
        beep
        disp('No model output found in this PRT.mat')
        return
        
    end
    
    % configure some variables
    CV       = PRT.model(modelid).input.cv_mat;     % CV matrix
    n_folds  = size(CV,2);                      % number of CV folds
    
    % parralel code?
    if def_par.allow
        try
            matlabpool(def_par.ncore)
        catch
            warning('Could not use pool of Matlab processes!')
        end
    end
    
    % targets
    t = PRT.model(modelid).input.targets;
    
    % load data files and configure ID matrix
    [Phi_all,ID,fid] = prt_getKernelModel(PRT,prt_dir,modelid);

    %get number of classes
    if strcmpi(PRT.model(modelid).input.type,'classification')
        nc=max(unique(t));
    else
        nc=[];
    end
    fdata.nc = nc;
    
    % Find chunks in the data (e.g. temporal correlated samples)
    % -------------------------------------------------------------------------
    
    ids = PRT.fs(fid).id_mat(PRT.model(modelid).input.samp_idx,:);
    i=1;
    samp_g=unique(ids(:,1));%number of groups
    for gid = 1: length(samp_g)
        
        samp_s=unique(ids(ids(:,1)==samp_g(gid),2)); %number of subjects for specific group
        
        for sid = 1: length(samp_s)
            
            samp_m=unique(ids(ids(:,1)==samp_g(gid) & ids(:,2)==samp_s(sid),3)); %number of modality for specific group & subject
            
            for mid = 1:length(samp_m)
                
                samp_c=unique(ids(ids(:,1)==samp_g(gid) & ids(:,2)==samp_s(sid) & ids(:,3)==samp_m(mid),4)); %number of conditions for specific group & subject & modality
                
                for cid = 1:length(samp_c)
                    
                    samp_b=unique(ids(ids(:,1)==samp_g(gid) & ids(:,2)==samp_s(sid) & ids(:,3)==samp_m(mid) & ids(:,4)==samp_c(cid),5));  %number of blocks for specific group & subject & modality & conditions
                    
                    for bid = 1:length(samp_b)
                        
                        rg = find((ids(:,1) == samp_g(gid)) & ...
                            (ids(:,2) == samp_s(sid)) & ...
                            (ids(:,3) == samp_m(mid)) & ...
                            (ids(:,4) == samp_c(cid)) & ...
                            (ids(:,5) == samp_b(bid)));
                        
                        chunks{i} = rg;
                        
                        i=i+1;
                    end
                end
            end
        end
    end
    
 
    
    % Initialize counts
    % -------------------------------------------------------------------------
    switch PRT.model(modelid).output.fold(1).type
        case 'classifier'
            n_class = length(PRT.model(modelid).output.fold(1).stats.c_acc);
            total_greater_c_acc = zeros(n_class,1);
            total_greater_b_acc = 0;
            
        case 'regression'
            total_greater_corr = 0;
            total_greater_mse = 0;
            total_greater_nmse = 0;
            total_greater_r2 = 0;
    end
    
    % Run model with permuted labels
    % -------------------------------------------------------------------------
    if ~isfield(PRT.model(modelid).output,'permutation') || ...
        (isfield(PRT.model(modelid).output,'permutation') && flag) %Back to empty to save other perm param
        PRT.model(modelid).output.permutation=struct('fold',[]);
    end
    for p=1:n_perm
        
        disp(sprintf('Permutation %d out of %d >>>>>>',p,n_perm));
        
        % permute
        chunkperm=randperm(length(chunks));
        CVperm = zeros(size(CV));
        t_perm = zeros(length(t),1);
        for i=1:length(chunks)
            t_perm(chunks{i},1)= unique(PRT.model(modelid).input.targets(chunks{chunkperm(i)})); 
            CVperm(chunks{i},:) = CV(chunks{chunkperm(i)},:);
        end
                
        for f = 1:n_folds
            % configure data structure for prt_cv_fold
            fdata.ID      = ID;
            fdata.mid     = modelid;
            fdata.CV      = CVperm(:,f);
            fdata.Phi_all = Phi_all;
            fdata.t       = t_perm;
            
            % Nested CV for hyper-parameter optimisation or feature selection
            if isfield(PRT.model(modelid).input,'use_nested_cv')
                if PRT.model(modelid).input.use_nested_cv
                    [out] = prt_nested_cv(PRT, fdata);
                    PRT.model(modelid).input.machine.args = out.opt_param;
                end
            end
            
            [temp_model, targets] = prt_cv_fold(PRT,fdata);
            
            % save the weights per fold to further compute ranking distance
            if flag
                PRT.model(modelid).output.permutation(p).fold(f).alpha=temp_model.alpha;
                PRT.model(modelid).output.permutation(p).fold(f).pred=temp_model.predictions;
            end
            
            model.output.fold(f).predictions = temp_model.predictions;
            model.output.fold(f).targets     = targets.test;
            
        end
        
        % Model level statistics (across folds)
        t             = vertcat(model.output.fold(:).targets);
        m.type        = PRT.model(modelid).output.fold(1).type;
        m.predictions = vertcat(model.output.fold(:).predictions);
        perm_stats         = prt_stats(m,t,t);
        
        
        switch PRT.model(modelid).output.fold(1).type
            
            case 'classifier'
                
                permutation.b_acc(p)=perm_stats.b_acc;
                n_class = length(PRT.model(modelid).output.fold(1).stats.c_acc);
                
                if (perm_stats.b_acc >= PRT.model(modelid).output.stats.b_acc)
                    total_greater_b_acc=total_greater_b_acc+1;
                end
                
                for c=1:n_class
                    permutation.c_acc(c,p)=perm_stats.c_acc(c);
                    if (perm_stats.c_acc(c) >= PRT.model(modelid).output.stats.c_acc(c))
                        total_greater_c_acc(c)=total_greater_c_acc(c)+1;
                    end
                end
                
            case 'regression'
                permutation.corr(p)=perm_stats.corr;
                if (perm_stats.corr >= PRT.model(modelid).output.stats.corr)
                    total_greater_corr=total_greater_corr+1;
                end
                permutation.mse(p)=perm_stats.mse;
                if (perm_stats.mse <= PRT.model(modelid).output.stats.mse)
                    total_greater_mse=total_greater_mse+1;
                end
                permutation.nmse(p)=perm_stats.nmse;
                if (perm_stats.nmse <= PRT.model(modelid).output.stats.nmse)
                    total_greater_nmse=total_greater_nmse+1;
                end
                permutation.r2(p)=perm_stats.r2;
                if (perm_stats.r2 >= PRT.model(modelid).output.stats.r2)
                    total_greater_r2=total_greater_r2+1;
                end
                
                
        end
    end
    
    switch PRT.model(modelid).output.fold(1).type
        case 'classifier'
            
            pval_b_acc = total_greater_b_acc / n_perm;
            if pval_b_acc == 0
                pval_b_acc = 1./n_perm;
            end
            
            pval_c_acc=zeros(n_class,1);
            for c=1:n_class
                pval_c_acc(c) = total_greater_c_acc(c) / n_perm;
                if pval_c_acc(c) == 0
                    pval_c_acc(c) = 1./n_perm;
                end
            end
            
            permutation.pvalue_b_acc = pval_b_acc;
            permutation.pvalue_c_acc = pval_c_acc;
            
        case 'regression'
            
            pval_corr = total_greater_corr / n_perm;
            if pval_corr == 0
                pval_corr = 1./n_perm;
            end
            
            pval_mse = total_greater_mse / n_perm;
            if pval_mse == 0
                pval_mse = 1./n_perm;
            end
            
            pval_nmse = total_greater_nmse / n_perm;
            if pval_nmse == 0
                pval_nmse = 1./n_perm;
            end
            
            pval_r2 = total_greater_r2 / n_perm;
            if pval_r2 == 0
                pval_r2 = 1./n_perm;
            end
            
            permutation.pval_corr = pval_corr;
            permutation.pval_mse = pval_mse;
            permutation.pval_nmse = pval_nmse;
            permutation.pval_r2 = pval_r2;
    end
    
    
    
    %update PRT
    PRT.model(modelid).output.stats.permutation = permutation;
    
    % Save PRT containing machine output
    % -------------------------------------------------------------------------
    outfile = fullfile(path,'PRT.mat');
    disp('Updating PRT.mat.......>>')
    if spm_check_version('MATLAB','7') < 0
        save(outfile,'-V6','PRT');
    else
        save(outfile,'PRT');
    end
    disp('Permutation test done.')
end

end

