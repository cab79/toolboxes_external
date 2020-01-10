function output = prt_machine_sMKL_cla(d,args)
% Run L1-norm MKL - wrapper for simpleMKL
% FORMAT output = prt_machine_sMKL_cla(d,args)
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
options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
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
options.nbitermax=def.model.l1MKLmaxitr;  % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this
options.seuilitermax=10;           % value, for iterations lower than this one

options.miniter=0;                 % minimal number of iterations
options.verbosesvm=0;              % verbosity of inner svm algorithm
options.efficientkernel=0;         % use efficient storage of kernels

%------------------------------------------------------
% Sanity check
%------------------------------------------------------
SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a string (empty or otherwise)
    if ~isnumeric(args)
        error('prt_machine_sMKL_cla:MKLargsNotNumber',['Error: L1_MKL'...
            ' args should be a number. ' ...
            ' SOLUTION: Please change range']);
    end
    
    % check it is indeed a two-class classification problem
    uTL=unique(d.tr_targets(:));
    nC=numel(uTL);
    if nC>2
        error('prt_machine_sMKL_cla:problemNotBinary',['Error:'...
            ' This machine is only for two-class problems but the' ...
            ' current problem has ' num2str(nC) ' ! ' ...
            'SOLUTION: Please select another machine than ' ...
            'prt_machine_sMKL_cla']);
    end
    % check it is indeed labelled correctly (probably should be done 
    if ~all(uTL==[1 2]')
        error('prt_machine_sMKL_cla:LabellingIncorect',['Error:'...
            ' This machine needs labels to be in {1,2} ' ...
            ' but they are ' mat2str(uTL) ' ! ' ...
            'SOLUTION: Please relabel your classes by changing the '...
            ' ''tr_targets'' argument to prt_machine_sMKL_cla']);
    end
end

%--------------------------------------------------------------------------
% Run simpleMKL
%--------------------------------------------------------------------------
C_opt = args;
options.sigmainit = 1/size(d.train,2)*ones(1,size(d.train,2)); %initialize kernel weights

% change targets from 1/2 to -1/1 
tr_targets = d.tr_targets;
c1PredIdx  = tr_targets  ==1; 
tr_targets  (c1PredIdx)  = 1; %positive values = 1 
tr_targets  (~c1PredIdx) = -1; %negative values = 2 

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

[beta,alpha_sv,b,pos] = mklsvm(ktrain,tr_targets,C_opt,options,verbose);

alpha = zeros(length(d.tr_targets),1);
alpha(pos) = alpha_sv;

ktest_final = zeros(length(d.te_targets),length(d.tr_targets));

for i = 1:size(d.train,2)
    ktest_final = ktest_final + beta(i)*ktest(:,:,i);
end

func_val = (ktest_final*alpha)+b;

predictions = sign(func_val);


% Outputs
%--------------------------------------------------------------------------
% change predictions from 1/-1 to 1/2 
c1PredIdx               = predictions==1; 
predictions(c1PredIdx)  = 1; %positive values = 1 
predictions(~c1PredIdx) = 2; %negative values = 2 

output.predictions = predictions;
output.func_val    = func_val;
output.type        = 'classifier';
output.alpha       = alpha;
output.b           = b;
output.totalSV     = length(alpha_sv);
output.beta        = beta; %kernel weights

end
