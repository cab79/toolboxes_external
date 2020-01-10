function output = prt_machine_svm_bin(d,args)
% Run binary SVM - wrapper for libSVM
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

% Written by M.J.Rosa, J.Mourao-Miranda and J.Richiardi
% $Id$

SANITYCHECK=true; % can turn off for "speed". Expert only.

%Turn the value of the C hyper-parameter into the arguments format for LIBSVM
if ~ischar(args)
    def = prt_get_defaults('model');
    args = [def.libsvmargs, num2str(args)];
end

if SANITYCHECK==true
    % args should be a string (empty or otherwise)
    if ~ischar(args)
        error('prt_machine_svm_bin:libSVMargsNotString',['Error: libSVM'...
            ' args should be a string. ' ...
            ' SOLUTION: Please do XXX']);
    end
    
    % check we can reach the binary library
    if ~exist('svmtrain','file')
        error('prt_machine_svm_bin:libNotFound',['Error:'...
            ' libSVM svmtrain function could not be found !' ...
            ' SOLUTION: Please check your path.']);
    end
    % check it is indeed a two-class classification problem
    uTL=unique(d.tr_targets(:));
    nC=numel(uTL);
    if nC>2
        error('prt_machine_svm_bin:problemNotBinary',['Error:'...
            ' This machine is only for two-class problems but the' ...
            ' current problem has ' num2str(nC) ' ! ' ...
            'SOLUTION: Please select another machine than ' ...
            'prt_machine_svm_bin in XXX']);
    end
    % check it is indeed labelled correctly (probably should be done 
    if ~all(uTL==[1 2]')
        error('prt_machine_svm_bin:LabellingIncorect',['Error:'...
            ' This machine needs labels to be in {1,2} ' ...
            ' but they are ' mat2str(uTL) ' ! ' ...
            'SOLUTION: Please relabel your classes by changing the '...
            ' ''tr_targets'' argument to prt_machine_svm_bin']);
    end
    
    % check we are using the C-SVC (exclude types -s 1,2,3,4)
    if ~isempty(regexp(args,'-s\s+[1234]','once'))
        error('prt_machine_svm_bin:argsProblem:onlyCSVCsupport',['Error:'...
            ' This machine only supports a C-SVC formulation ' ...
            ' (''-s 0'' in the ''args'' parameter), but the args ' ...
            ' supplied are ''' args ''' ! ' ...
            'SOLUTION: Please change the offending part of args to '...
            '''-s 0''']);
    end
    
    % check we are using linear or precomputed kernels
    % (exclude types -t 1,2,3)
    if ~isempty(regexp(args,'-t\s+[123]','once'))
        error('prt_machine_svm_bin:argsProblem:onlyLinOrPrecomputeSupport',...
            ['Error: This machine only supports linear or precomputed ' ...
            'kernels (''-t 0/4'' in the ''args'' parameter), but the args ' ...
            ' supplied are ''' args ''' ! ' ...
            'SOLUTION: Please change the offending part of args to '...
            '''-t 0'' or ''-t 4'' as intended']);
    end
    
end


% Run SVM
%--------------------------------------------------------------------------
nlbs  = length(d.tr_targets);
allids_tr = (1:nlbs)';

model = svmtrain(d.tr_targets,[allids_tr d.train{:}],args);

% check if training succeeded:
if isempty(model)
    if (ischar(args))
        args_str = args;
    else
        args_str = '';
    end
    error('prt_machine_svm_bin:libSVMsvmtrainUnsuccessful',['Error:'...
        ' libSVM svmtrain function did not run properly!' ...
        ' This could be a problem with the supplied function arguments'...
        ' ' args_str '']);
end


% Get SV coefficients (alpha) in the original order and the bias term (b) 
sgn   = -1*(2 * model.Label(1) - 3); %variable to account for label convention in PRoNTo
alpha = get_alpha(model,nlbs,sgn);
b     = -model.rho *sgn;

% compute prediction directly rather than using svmpredict, which does
% not allow empty test labels
if iscell(d.test)
    func_val = cell2mat(d.test)*alpha+b;
else
    func_val = d.test*alpha+b;
end

% compute hard decisions
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
output.totalSV     = model.totalSV;

end

% Get SV coefficients
%--------------------------------------------------------------------------
function alpha = get_alpha(model,n,sgn)
% needs a function because examples can be re-ordered by libsvm
alpha = zeros(n,1);

for i = 1:model.totalSV
    ind        = model.SVs(i);
    alpha(ind) = model.sv_coef(i);
end

alpha = sgn*alpha;

end

