function output = prt_machine_RT_bin(d,args)
% Run binary Ensemble of Regression Tree - wrapper for Pierre Geurt's
% RT code
% FORMAT output =  prt_machine_RT_bin(d,args)
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
%    args    - vector of RT arguments
%       args(1) - number of trees (default: 501)
% Output:
%    output  - output of machine (struct).
%     * Mandatory fields:
%      .predictions - predictions of classification or regression [Nte x D]
%     * Optional fields:
%      .func_val - value of the decision function
%      .type     - which type of machine this is (here, 'classifier')
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

%--------------------------------------------------------------------------
% Written by J.Richiardi
% $Id$

% FIXME: support for multiple kernels / feature representations
% is not yet tested, there might be transposition or dimensionality errors.

% TODO: this machine supports regression AND classification
% TODO: this machine supports multi-class

SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a vector (empty or otherwise)
    if ~isvector(args)
        error('prt_machine_RT_bin:RTargsNotVec',['Error: RT'...
            ' args should be a vector. ' ...
            ' SOLUTION: Please check your code. ']);
    end
    
    % check we can reach the binary library
    if ~exist('rtenslearn_c','file')
        error('prt_machine_RT_bin:libNotFound',['Error:'...
            ' RT function rtenslearn_c could not be found !' ...
            ' SOLUTION: Please check your path.']);
    end
    
    % check it is indeed a two-class classification problem
    uTL=unique(d.tr_targets(:)); % unique training labels
    nC=numel(uTL);
    if nC>2
        error('prt_machine_RT_bin:problemNotBinary',['Error:'...
            ' This machine is only for two-class problems but the' ...
            ' current problem has ' num2str(nC) ' ! ' ...
            'SOLUTION: Please select another machine than ' ...
            'prt_machine_RT_bin in XXX']);
    end
    
    % check it is indeed labelled correctly (probably should be done
    % above?)
    if ~all(uTL==[1 2]')
        error('prt_machine_RT_bin:LabellingIncorect',['Error:'...
            ' This machine needs labels to be in {1,2} ' ...
            ' but they are ' mat2str(uTL) ' ! ' ...
            'SOLUTION: Please relabel your classes by changing the '...
            ' ''tr_lbs'' argument to prt_machine_RT_bin']);
    end
    
    % check we are not setting a ridiculous number of trees
    if isempty(args)
        args(1)=501;
        disp('prt_machine_RT_bin: defaulting to 501 trees');
    else
        if args(1)>10000;
            error('prt_machine_RT_bin:argsProblem:maxTrees',['Error:'...
                ' Setting a high number of trees is not supported ' ...
                ' without some modifications of the wrapper code. ' ...
                ' Expert only! ' ...
                'SOLUTION: Please change the offending args to '...
                'a value less than 10000.']);
        end
    end
end % SANITYCHECK


% Run RT
%--------------------------------------------------------------------------
rtParams=init_rf(); % random forests
rtParams.nbterms=args(1); % number of trees
tridx=int32(1:numel(d.tr_targets));  % (WARNING: int32 format is mandatory)
verbose=1;   % TODO: make this a machine arg

[output.func_val output.w trees]=rtenslearn_c(single(d.train{1}),...
    single(d.tr_targets),tridx,[],rtParams,single(d.test{1}),verbose);

% check if training succeeded:
if isempty(output)
    error('prt_machine_RT_bin:RTtrainUnsuccessful',['Error:'...
        ' RT rtenslearn_c function did not run properly!' ...
        ' This could be a problem with the supplied function arguments'...
        mat2str(args)]);
end

% compute hard decisions
output.predictions=round(output.func_val);

if d.use_kernel==false
    % normalise importance to norm 1
    output.w=output.w/norm(output.w,1);
else
    % do nothing - we can't compute primal weights from inside here
end

% prepare output
output.type        = 'classifier';

end
