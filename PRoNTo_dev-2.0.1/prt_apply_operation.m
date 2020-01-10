function out = prt_apply_operation(PRT, in, opid)
% function to apply a data operation to the training, test and 
% in.train:      training data
% in.tr_id:      id matrix for training data
% in.use_kernel: are the data in kernelised form
% in.tr_targets: training targets (optional field)
% in.pred_type:  'classification' or 'regression' (required for tr_targets)
%
% A test set may also be specified, which require the following fields:
% in.test:       test data
% in.testcov:    test covariance (only if use_kernel = true)
% in.te_targets: test targets
% in.te_id:      id matrix for test data
%
% opid specifies the operation to apply, where:
%    1 = Temporal Compression
%    2 = Sample averaging (average samples for each subject/condition)
%    3 = Mean centre features over subjects
%    4 = Divide data vectors by their norm
%    5 = Perform a GLM (fMRI only)
%
% N.B: - all operations are applied independently to training and test
%        partitions
%      - see Chu et. al (2011) for mathematical descriptions of operations
%        1 and 2 and Shawe-Taylor and Cristianini (2004) for a description
%        of operation 3.
%
% References:
% Chu, C et al. (2011) Utilizing temporal information in fMRI decoding: 
% classifier using kernel regression methods. Neuroimage. 58(2):560-71.
% Shawe-Taylor, J. and Cristianini, N. (2004). Kernel methods for Pattern
% analysis. Cambridge University Press.
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

% copy input fields to output
out = in;

for d = 1:length(in.train)
    switch opid
        case 1 
            % temporal compression
            % --------------------
            % Training data
            Ptr = compute_tc_mat(in.tr_id);
            if in.use_kernel
                out.train{d} = Ptr*in.train{d}*Ptr';
            else
                out.train{d} = Ptr*in.train{d};
            end
            out.tr_id = round(Ptr*in.tr_id);
            if isfield(in,'tr_targets')
                out.tr_targets = Ptr*in.tr_targets;
                if strcmpi(in.pred_type,'classification');
                    out.tr_targets = round(out.tr_targets);
                end
            end
            
            % Test data 
            if isfield(in,'test')
                Pte = compute_tc_mat(in.te_id);
                if in.use_kernel
                    out.test{d}     = Pte*in.test{d}*Ptr';
                    out.testcov{d}  = Pte*in.testcov{d}*Pte';
                else
                    out.test{d}  = Pte*in.test{d};
                end
                out.te_id = round(Pte*in.te_id);
                if isfield(in,'te_targets')
                    out.te_targets = Pte*in.te_targets;
                    if strcmpi(in.pred_type,'classification');
                        out.te_targets = round(out.te_targets);
                    end
                end
            end
            
        case 2  
            % sample averaging
            % ----------------
            % Training data
            Ptr = compute_sa_mat(in.tr_id);
            if in.use_kernel
                out.train{d} = Ptr*in.train{d}*Ptr';
            else
                out.train{d} = Ptr*in.train{d};
            end
            out.tr_id = round(Ptr*in.tr_id);
            if isfield(in,'tr_targets')
                out.tr_targets = Ptr*in.tr_targets;
                if strcmpi(in.pred_type,'classification');
                    out.tr_targets = round(out.tr_targets);
                end
            end
            
            % Test data
            if isfield(in,'test')
                Pte = compute_sa_mat(in.te_id);
                if in.use_kernel
                    out.test{d}     = Pte*in.test{d}*Ptr';
                    out.testcov{d}  = Pte*in.testcov{d}*Pte';
                else
                    out.test{d}  = Pte*in.test{d};
                end
                out.te_id = round(Pte*in.te_id);
                if isfield(in,'te_targets')
                    out.te_targets = Pte*in.te_targets;
                    if strcmpi(in.pred_type,'classification');
                        out.te_targets = round(out.te_targets);
                    end
                end
            end
            
        case 3 
            % mean centre features over subjects
            % ----------------------------------
            if ~isfield(in,'test') 
                % No test data
                if in.use_kernel
                    out.train{d} = prt_centre_kernel(in.train{d});
                else
                    m = mean(in.train{d});
                    %out.train{d} = in.train{d} - repmat(m,size(in.train{d},2),1);
                    out.train{d} = zeros(size(in.train{d}));
                    for r = 1:size(in.train{d},1)
                        out.train{d}(r,:) = in.train{d}(r,:) - m;
                    end
                end
            else % Test data supplied
                 if in.use_kernel
                    [out.train{d}, out.test{d}, out.testcov{d}] = ...
                        prt_centre_kernel(in.train{d},in.test{d},in.testcov{d});
                else
                    m = mean(in.train{d});
                    %out.train{d} = in.train{d} - repmat(m,size(in.train{d},2),1);
                    %out.test{d}  = in.test{d} - repmat(m,size(in.test{d},2),1);
                    out.train{d} = zeros(size(in.train{d}));
                    for r = 1:size(in.train{d},1)
                        out.train{d}(r,:) = in.train{d}(r,:) - m;
                    end
                    out.test{d} = zeros(size(in.test{d}));
                    for r = 1:size(in.test{d},1)
                        out.test{d}(r,:)  = in.test{d}(r,:) - m;
                    end
                end
                out.te_id = in.te_id;
            end
            out.tr_id = in.tr_id;
            if isfield(in,'tr_targets')
                out.tr_targets = in.tr_targets;
            end
            if isfield(in,'te_targets')
                out.te_targets = in.te_targets;
            end
            
        case 4 
            % divide each feature vector by its norm
            % --------------------------------------
            % in this case, the operation is applied independently to
            % each data vector, so it is safe (and convenient) to apply
            % the operation to the whole kernel at once
            if ~isfield(in,'test')
                % No test data
                if in.use_kernel
                    Phi = prt_normalise_kernel(in.train{d});
                    tr = 1:size(in.train{d},1);
                    out.train{d} = Phi(tr,tr);
                else
                    out.train{d} = zeros(size(in.train{d}));
                    for r = 1:size(in.train{d})
                        out.train{d}(r,:) = in.train{d}(r,:) / norm(in.train{d}(r,:));
                    end
                end
            else % Test data
                
                if in.use_kernel
                    Phi = [in.train{d}, in.test{d}'; in.test{d}, in.testcov{d}];
                    Phi = prt_normalise_kernel(Phi);
                    
                    tr = 1:size(in.train{d},1);
                    te = (1:size(in.test{d},1))+max(tr);
                    out.train{d}    = Phi(tr,tr);
                    out.test{d}     = Phi(te,tr);
                    out.testcov{d}  = Phi(te,te);
                else
                    out.train{d} = zeros(size(in.train{d}));
                    for r = 1:size(in.train{d})
                        out.train{d}(r,:) = in.train{d}(r,:) / norm(in.train{d}(r,:));
                    end
                    out.train{d} = zeros(size(in.test{d}));
                    for r = 1:size(in.test{d})
                        out.test{d}(r,:) = in.test{d}(r,:) / norm(in.test{d}(r,:));
                    end
                end
                out.te_id = in.te_id;
            end
            out.tr_id = in.tr_id;
            if isfield(in,'tr_targets')
                out.tr_targets = in.tr_targets;
            end
            if isfield(in,'te_targets')
                out.te_targets = in.te_targets;
            end
            
        case 5 
            % perform a GLM
            % -------------
            if ~isfield(in,'tr_cov')
                error('prt_apply_operation:NoCovariates',...
                'No covariates found to perform requested GLM');
            end
            if ~isfield(in,'test') 
                % No test data
                if in.use_kernel
                    out.train{d} = prt_remove_confounds(in.train{d},...
                        [in.tr_cov,ones(size(in.tr_cov,1),1)]);
                else
                    error('prt_apply_operation:GLMnonKernel',...
                'GLM not implemented for non-kernel methods');
                end
            else % Test data: for now does NOT take train/test division
                Phi = [in.train{d}, in.test{d}'; in.test{d}, in.testcov{d}];
                 if in.use_kernel
                    C = [in.tr_cov;in.te_cov];
                    C = [C, ones(size(C,1),1)];
                    [Phi] = prt_remove_confounds(Phi,C);
                    tr = 1:size(in.train{d},1);
                    te = (1:size(in.test{d},1))+max(tr);
                    out.train{d}    = Phi(tr,tr);
                    out.test{d}     = Phi(te,tr);
                    out.testcov{d}  = Phi(te,te);
                else
                    error('prt_apply_operation:GLMnonKernel',...
                     'GLM not implemented for non-kernel methods');
                end
                out.te_id = in.te_id;
            end
            out.tr_id = in.tr_id;
            if isfield(in,'tr_targets')
                out.tr_targets = in.tr_targets;
            end
            if isfield(in,'te_targets')
                out.te_targets = in.te_targets;
            end
            
                        
        otherwise
            error('prt_apply_operation:UnknownOperationSpecified',...
                'Unknown operation requested');
    end
end

%out.use_kernel = in.use_kernel;
%if isfield(in,'pred_type');
%    out.pred_type  = in.pred_type;
%end
end

% -------------------------------------------------------------------------
% Private Functions
% -------------------------------------------------------------------------

function P = compute_tc_mat(ID)
% function to compute the block averaging matrix (P) necessary to apply
% temporal compression

% give each block a unique id
IDc = zeros(size(ID,1),1);
C = {}; 
ccount = 0; 
lastid = zeros(1,5);
for c = 1:size(ID,1)
    currid = ID(c,1:5);  
    if any(lastid ~= currid)
        ccount = ccount + 1;
    end
    lastid = currid;
    IDc(c) = ccount;
end

% Compute sample averaging matrix
cids  = unique(IDc);
cnums = histc(IDc,cids);
C = cell(length(cnums),1);
for c = 1:length(cnums)
    C{c} = 1/cnums(c) .* ones(1,cnums(c));
end
P = blkdiag(C{:});
end

function P = compute_sa_mat(ID)
% function to compute the block averaging matrix (P) necessary to apply
% temporal compression

% give each subject a unique id
IDs = zeros(size(ID,1),1);
ccount = 0; 
lastid = zeros(1,2);
for s = 1:size(ID,1)
    currid = ID(s,1:2);  
    if any(lastid ~= currid)
        ccount = ccount + 1;
    end
    lastid = currid;
    IDs(s) = ccount;
end

subs = unique(IDs);

P = [];
for s = 1:length(subs)
    sidx = IDs == subs(s);
    conds = unique(ID(sidx,4));
    for c = 1:length(conds)
        p = (IDs == s & ID(:,4) == conds(c))';
        P = [P; 1./sum(p) * double(p)];
    end
end
P = double(P);
end

