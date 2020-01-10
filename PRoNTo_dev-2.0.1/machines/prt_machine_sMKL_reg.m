function output = prt_machine_sMKL_reg(d,args)
% Run L1-norm MKL - wrapper for simpleMKL
% FORMAT output = prt_machine_sMKL_reg(d,args)
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
%    args     - simpleMKL arguments
% Output:
%    output  - output of machine (struct).
%     * Mandatory fields:
%      .predictions - predictions of classification or regression [Nte x D]
%     * Optional fields:
%      .func_val - value of the decision function
%      .type     - which type of machine this is (here, 'classifier')
%      .
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Mourao-Miranda 

def = prt_get_defaults;

%------------------------------------------------------
% configure simpleMKL options
%------------------------------------------------------
verbose=0;
options.algo='svmreg'; % Choice of algorithm in mklsvm can be either
% 'svmclass' or 'svmreg'

%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=0; % use variation of weights for stopping criterion
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion
options.stopdualitygap=1; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
% are set to zero
options.lambdareg = 1e-8;          % ridge added to kernel matrix

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base
% variable in the reduced gradient method
options.nbitermax=def.model.l1MKLmaxitr;;             % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this
options.seuilitermax=10;           % value, for iterations lower than this one

options.miniter=0;                 % minimal number of iterations
options.verbosesvm=0;              % verbosity of inner svm algorithm
options.efficientkernel=0;         % use efficient storage of kernels
options.svmreg_epsilon=0.01;

% Run simpleMKL
%--------------------------------------------------------------------------
C_opt = args;
options.sigmainit = 1/size(d.train,2)*ones(1,size(d.train,2)); %initialize kernel weights
m = mean(d.tr_targets);  % mean of the training data
tr_targets = d.tr_targets - m; % mean centre targets

%reshape previously normalized kernel
ktrain = zeros(size(d.train{1},1),size(d.train{1},1),size(d.train,2));
ktest = zeros(size(d.test{1},1),size(d.train{1},1),size(d.train,2));
for k = 1:size(d.train,2)
    if sum(sum(isnan(d.train{k})))==0;
        ktrain(:,:,k) =  d.train{k}  ;
    end
    if sum(sum(isnan(d.test{k}))) ==0
        ktest(:,:,k) =  d.test{k}   ;
    end
end

[beta,alpha_sv,b,pos,history,obj,status] = mklsvm(ktrain,tr_targets,C_opt,options,verbose);

alpha = zeros(length(d.tr_targets),1);
alpha(pos) = alpha_sv;

ktest_final = zeros(length(d.te_targets),length(d.tr_targets));

for i = 1:size(d.train,2)
    ktest_final = ktest_final + beta(i)*ktest(:,:,i);
end

func_val = ((ktest_final*alpha)+b)+m; % add mean from the training set

predictions = func_val;

% Outputs
%-------------------------------------------------------------------------
output.predictions = predictions;
output.func_val    = func_val;
output.type        = 'regression';
output.alpha       = alpha;
output.b           = b;
output.totalSV     = length(alpha_sv);
output.beta        = beta; %kernel weights

end