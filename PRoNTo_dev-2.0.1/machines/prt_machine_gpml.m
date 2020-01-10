function output = prt_machine_gpml(d,args)
% Run Gaussian process model - wrapper for gpml toolbox
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

% configure default parameters for GP optimisation
meanfunc  = @meanConstcell;
covfunc   = @covLINkcell; 
maxeval   = -20;
if strcmp(d.pred_type,'classification')
    mode = 'classifier'; % it's good to be consistent!
    likfunc   = @likErf;
    inffunc   = @prt_infEP;
else
    mode = 'regression';
    likfunc   = @likGauss;
    inffunc   = @prt_infExact;
    mtr       = mean(d.tr_targets);      % mean of the training data
end

% Error checks
% -------------------------------------------------------------------------
SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a string (empty or otherwise)
    if ~ischar(args)
        error('prt_machine_gpml:argsNotString',['Error: gpml'...
            ' args should be a string. ' ...
            ' SOLUTION: Please do XXX']);
    end
    
    % check we can reach the binary library
    if ~exist('prt_gp','file')
        error('prt_machine_gpml:libNotFound',['Error:'...
            ' ''prt_gp'' function could not be found !' ...
            ' SOLUTION: Please check your path.']);
    end
    % check whether it is a two-class classification problem
    uTL=unique(d.tr_targets(:));
    k=numel(uTL); % number of classes
    if strcmp(mode,'classifier') && k > 2
        warning('prt_machine_gpml:classificationWithMoreThanTwoClasses',...
               ['Classification specified with > 2 classes. ',...
                'Defaulting to multiclass Laplace approximation.']);
        output = prt_machine_gpclap(d,args);
        return;
    end
    % are we using a kernel ? 
    if ~d.use_kernel
        error('prt_machine_gpml:useKernelIsFalse',['Error:'...
            ' This machine is currently only implemented for kernel data ' ...
            'SOLUTION: Please set use_kernel to true']);
    end
end

% parse input arguments
% -------------------------------------------------------------------------
% hyperparameters
if ~isempty(regexp(args,'-h','once'))
    optimise_theta = true;
    eargs = regexp(args,'-f\s+[0-9]*','match');
    if ~isempty(eargs)
        eargs = regexp(cell2mat(eargs),'-f\s+','split');
        maxeval  = str2num(['-',cell2mat(eargs(2))]);
    end
else
    optimise_theta = false;
end
% likelihood function
largs = regexp(args,'-l\s+[a-zA-Z0-9_]*','match');
if ~isempty(largs)
    largs = regexp(cell2mat(largs),'-l\s+','split');
    likfunc = str2func(cell2mat(largs(2)));
    if strcmpi(cell2mat(largs(2)),'Erf')
        likfunc   = @likErf;
    end
end
% covariance function
cargs = regexp(args,'-c\s+[a-zA-Z0-9_]*','match');
if ~isempty(cargs)
    cargs = regexp(cell2mat(cargs),'-c\s+','split');
    covfunc = str2func(cell2mat(cargs(2)));
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
    inffunc = str2func(cell2mat(iargs(2)));
end
% priors
if ~isempty(regexp(args,'-p','once'))
    disp('Empirical priors specified. Using MAP for hyperparameters')
    priors = prt_gp_priors(func2str(covfunc));
    map = true;
else
    map = false;
end

% Set default hyperparameters
% -------------------------------------------------------------------------
nhyp = str2num([feval(covfunc); feval(likfunc); feval(meanfunc)]);
if nhyp(1) > 0
    hyp.cov = zeros(nhyp(1),1);
end
if nhyp(2) > 0
    hyp.lik = zeros(nhyp(2),1);
end
if nhyp(3) > 0 
    hyp.mean = zeros(nhyp(3),1);
end

% Assemble data matrices
% -------------------------------------------------------------------------
% handle the glm as a special case (for now)
if strcmpi(func2str(covfunc),'covLINglm') || strcmpi(func2str(covfunc),'covLINglm_2class')
    % configure covariances
    K   = [d.train(:)'   {d.tr_param}];
    Ks  = [d.test(:)'    {d.te_param}];
    Kss = [d.testcov(:)' {d.te_param}];
    
    % get default hyperparamter values
    hyp.cov = log(prt_glm_design);
    
    [tmp1 tmp2 tmp3 tr_lbs] = prt_glm_design(hyp.cov, d.tr_param);
    [tmp1 tmp2 tmp3 te_lbs] = prt_glm_design(hyp.cov, d.te_param);   
else
    % configure covariances
    K   = d.train;
    Ks  = d.test;
    Kss = d.testcov;
        
    tr_lbs = d.tr_targets;
    te_lbs = d.te_targets;
end

% configure targets
if strcmp(mode,'classifier')
    % convert targets to +1/-1
    y = -1*(2 * tr_lbs - 3);
else
    y = tr_lbs - mtr;
end
    
% Train and test GP model
% -------------------------------------------------------------------------
% train
if optimise_theta
    if map
        [hyp,nlmls] = minimize(hyp, @prt_gp_map, maxeval, inffunc, meanfunc, covfunc, likfunc, K, y, priors);
    else
        [hyp nlmls] = minimize(hyp, @prt_gp, maxeval, inffunc, meanfunc, covfunc, likfunc, K, y);
    end
else
    nlmls = prt_gp(hyp, inffunc, meanfunc, covfunc, likfunc, K, y);
end

% make predictions
[ymu ys2 fmu fs2 lp post] = prt_gp(hyp, inffunc, meanfunc, covfunc, likfunc,K, y, Ks, zeros(size(Ks{1},1),1), Kss);

% Outputs
% -------------------------------------------------------------------------
if strcmp(mode,'classifier')
    p = exp(lp);
    output.predictions = (1-real(p > 0.5)) + 1;
    output.func_val    = p;
else % regression
    output.predictions = ymu + mtr;
    output.func_val    = output.predictions;
end
output.type        = mode;
output.loghyper    = hyp;
output.mu          = ymu;
output.s2          = ys2;
output.nlml        = min(nlmls);
output.tr_targets  = tr_lbs;
output.te_targets  = te_lbs;
output.alpha       = post.alpha;
%output.sW          = post.sW;
%output.L           = post.L;

end

