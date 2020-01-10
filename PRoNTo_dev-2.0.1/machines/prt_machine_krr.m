function output = prt_machine_krr(d,args)
% Kernel ridge regression
% FORMAT output = prt_machine_svm_bin(d,args)
% Inputs:
%   d         - structure with data information, with mandatory fields:
%     .train      - training data (cell array of matrices of row vectors,
%                   each [Ntr x D]). each matrix contains one representation
%                   of the data. This is useful for approaches such as
%                   multiple kernel learning.
%     .test       - testing data  (cell array of matrices row vectors, each
%                   [Nte x D])
%     .tr_targets - training labels (for classification) or values (for
%                   regression) (column vector, [Ntr x 1])
%     .use_kernel - flag, is data in form of kernel matrices (true) of in
%                form of features (false)
%    args     - libSVM arguments
% Output:
%    output  - output of machine (struct).
%     * Mandatory fields:
%      .predictions - predictions of classification or regression [Nte x D]
%     * Optional fields:
%      .func_val - value of the decision function
%      .type     - which type of machine this is (here, 'classifier')
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Carlton Chu
% $Id$


output.type=d.pred_type;

SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a string (empty or otherwise)
    if ~isnumeric(args)
        error('prt_machine_krr:krrRegNotNumeric',['Error: krr'...
            ' regularization should be a number. ' ...
            ' SOLUTION: Please use a number']);
    end
    K=d.train{1};
    [n m]=size(K);
    if n~=m
        error('prt_machine_krr:kernelSize',['Error: krr'...
            ' training kernel should be square ' ...
            ' SOLUTION: do the right thing']);
    end
    
    % Run KRR
    %----------------------------------------------------------------------
    m = mean(d.tr_targets);                     % mean of the training data
    t = d.tr_targets - m;                             % mean centre targets
    w = prt_KRR(K,t,args);
       
    output.predictions=d.test{1}*w + m;    % add mean from the training set 
    output.func_val=output.predictions;
    output.alpha=w;
  
end

end
