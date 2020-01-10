function output = prt_machine_gpr(d,args)
% Run Gaussian process regression - meta-wrapper for regression with gpml  
% FORMAT output = prt_machine_gpml(d,args)
% Inputs:
%   d         - structure with data information, with mandatory fields:
%     .train      - training data (cell array of matrices of row vectors,
%                   each [Ntr x D]). each matrix contains one representation
%                   of the data. This is useful for approaches such as
%                   multiple kernel learning.
%     .test       - testing data  (cell array of matrices row vectors, each
%                   [Nte x D])
%     .testcov    - testing covariance (cell array of matrices row vectors,
%                   each [Nte x Nte])
%     .tr_targets - training labels (for classification) or values (for
%                   regression) (column vector, [Ntr x 1])
%     .use_kernel - flag, is data in form of kernel matrices (true) or in 
%                   form of features (false)
%    args     - argument string, where
%       -h         - optimise hyperparameters (otherwise don't)
%       -f iter    - max # iterations for optimiser (ignored if -h not set)
%       -l likfun  - likelihood function:
%                       'likErf' - erf/probit likelihood (binary only)
%       -c covfun  - covariance function:
%                       'covLINkcell' - simple dot product
%                       'covLINglm'   - construct a GLM
%       -m meanfun - mean function:
%                       'meanConstcell' - suitable for dot product
%                       'meanConstglm'  - suitable for GLM
%       -i inffun  - inference function:
%                       'prt_infEP' - Expectation Propagation
%    experimental args (use at your own risk):
%       -p         - use priors for the hyperparameters. If specified, this
%                    indicates that a maximum a posteriori (MAP) approach
%                    will be used to set covariance function
%                    hyperparameters. The priors are obtained by calling
%                    prt_gp_priors('covFuncName')
%
%       N.B.: for the arguments specifying functions, pass in a string, not
%       a function handle. This script will generate a function handle
% 
% Output:
%    output  - output of machine (struct).
%     * Mandatory fields:
%      .predictions - predictions of classification or regression [Nte x D]
%     * Optional fields:
%      .type     - which type of machine this is (here, 'classifier')
%      .func_val - predictive probabilties
%      .mu       - test latent means
%      .s2       - test latent variances
%      .loghyper - log hyperparameters
%      .nlml     - negative log marginal likelihood
%      .alpha    - GP weighting coefficients
%      .sW       - likelihood matrix (see Rasmussen & Williams, 2006)
%      .L        - Cholesky factor
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

% set default gp paramters (N.B.: these are strings, not function handles)
meanfunc  = 'meanConstcell';
covfunc   = 'covLINkcell'; 
maxeval   = '100';
likfunc   = 'likGauss';
inffunc   = 'prt_infExact';

% parse input arguments (i.e. check for non-default options)
% -------------------------------------------------------------------------
% hyperparameters
if ~isempty(regexp(args,'-h','once'))
    opt = ' -h ';
    eargs = regexp(args,'-f\s+[0-9]*','match');
    if ~isempty(eargs)
        eargs = regexp(cell2mat(eargs),'-f\s+','split');
        maxeval  = [cell2mat(eargs(2))];
    end
else
    opt = '';
end
% likelihood function
largs = regexp(args,'-l\s+[a-zA-Z0-9_]*','match');
if ~isempty(largs)
    largs = regexp(cell2mat(largs),'-l\s+','split');
    likfunc = str2func(cell2mat(largs(2)));
    if strcmpi(cell2mat(largs(2)),'Gauss')
        likfunc  = 'likGauss';
    end
end
% covariance function
cargs = regexp(args,'-c\s+[a-zA-Z0-9_]*','match');
if ~isempty(cargs)
    cargs = regexp(cell2mat(cargs),'-c\s+','split');
    covfunc = cell2mat(cargs(2));
end
% mean function
margs = regexp(args,'-m\s+[a-zA-Z0-9_]*','match');
if ~isempty(margs)
    margs = regexp(cell2mat(margs),'-m\s+','split');
    meanfunc = str2func(cell2mat(margs(2)));
end
% inference function
iargs = regexp(args,'-i\s+[a-zA-Z0-9_]*','match');
if ~isempty(iargs)
    iargs = regexp(cell2mat(iargs),'-i\s+','split');
    inffunc = cell2mat(iargs(2));
end
% priors
if ~isempty(regexp(args,'-p','once'))
    map = ' -p ';
else
    map = '';
end

% construct argument string for prt_machine_gpml
args = ['-l ',likfunc,' -c ',covfunc,' -m ',meanfunc,...
       ' -i ',inffunc,' ',opt,' ','-f ',maxeval,map];

% do the regression
output = prt_machine_gpml(d, args);
end

