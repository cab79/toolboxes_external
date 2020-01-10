function output = prt_machine_gpclap(d,args)
% Run multiclass Gaussian process classification (Laplace approximation)
% FORMAT output = prt_machine_gpclap(d,args)
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
%       -c covfun  - covariance function:
%                       'covLINkcell' - simple dot product
%                       'covLINglm'   - construct a GLM
%    experimental args (use at your own risk):
%       -p         - use priors for the hyperparameters. If specified, this
%                    indicates that a maximum a posteriori (MAP) approach
%                    will be used to set covariance function
%                    hyperparameters. The priors are obtained 
%                    by calling prt_gp_priors('covFuncName')
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
%      .loghyper - log hyperparameters
%      .nlml     - negative log marginal likelihood
%      .mu       - test latent means
%      .s2       - test latent variances
%      .alpha    - GP weighting coefficients
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J Ashburner and A Marquand
% $Id$

% Error checks
% -------------------------------------------------------------------------
SANITYCHECK=true; % can turn off for "speed". Expert only.

if SANITYCHECK==true
    % args should be a string (empty or otherwise)
    if ~ischar(args)
        error('prt_machine_gpclap:libSVMargsNotString',['Error: gpml'...
            ' args should be a string. ' ...
            ' SOLUTION: Please do XXX']);
    end
    % are we using a kernel ?
    if ~d.use_kernel
        error('prt_machine_gpclap:useKernelIsFalse',['Error:'...
            ' This machine is currently only implemented for kernel data ' ...
            'SOLUTION: Please set use_kernel to true']);
    end
end

% configure default parameters for GP optimisation
covfunc   = @covLINkcell;
mode      = 'classifier';

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
% covariance function
cargs = regexp(args,'-c\s+[a-zA-Z0-9_]*','match');
if ~isempty(cargs)
    cargs = regexp(cell2mat(cargs),'-c\s+','split');
    covfunc = str2func(cell2mat(cargs(2)));
end
% priors
if ~isempty(regexp(args,'-p','once'))
    disp('Empirical priors specified. Using MAP for hyperparameters')
    priors = prt_gp_priors(func2str(covfunc));
    map = true;
else
    map = false;
end

% Set default hyperparameters and objective function
% -------------------------------------------------------------------------
nhyp = str2num(feval(covfunc));
if nhyp(1) > 0
    hyp = zeros(nhyp(1),1);
end
if map
    objfunc = @gp_objfun_map;
else
    objfunc = @gp_objfun;
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
    hyp = log(prt_glm_design);
    
    [tmp1 tmp2 tmp3 tr_lbs] = prt_glm_design(hyp, d.tr_param);
    [tmp1 tmp2 tmp3 te_lbs] = prt_glm_design(hyp, d.te_param);    
else
    % configure covariances
    K   = d.train;
    Ks  = d.test;
    Kss = d.testcov;
    
    tr_lbs = d.tr_targets;
    te_lbs = d.te_targets;
end 

% create one-of-k labels
k = max(unique(tr_lbs));
n = length(tr_lbs);
Y = zeros(n,k);
for j = 1:n 
    Y(j,tr_lbs(j)) = 1;
end

% Train and test GP model
% -------------------------------------------------------------------------
% train
if optimise_theta
    nh = numel(hyp);
    if map
        objfunc = @gp_objfun_map;
    end
    hyp = spm_powell(hyp,eye(nh),ones(nh,1)*0.05,objfunc,Y,K,covfunc); 
end
% compute marginal likelihood and posterior parameters
[f lml alpha] = gp_lap_multiclass(hyp,covfunc,K,Y);

% make predictions
[p mu sigma]  = gp_pred_lap_multiclass(hyp,K,Y,covfunc,Ks,Kss);
[maxp pred]   = max(p,[],2);

% Outputs
% -------------------------------------------------------------------------
output.predictions = pred;
output.type        = mode;
output.func_val    = p;
output.tr_targets  = tr_lbs;
output.te_targets  = te_lbs; 
output.mu          = mu;
output.sigma       = sigma;
output.loghyper    = hyp;
output.nlml        = -lml;
output.alpha       = alpha;
end

% -------------------------------------------------------------------------
% Private functions
% -------------------------------------------------------------------------

function E = gp_objfun(logtheta,t,X,covfunc)
% Objective function to minimise

[f,F]   = gp_lap_multiclass(logtheta,covfunc,X,t);
E = -F; %+ 1e-6*(logtheta'*logtheta);
end

% -------------------------------------------------------------------------
function E = gp_objfun_map(logtheta,t,X,covfunc)
% Objective function to minimise in a MAP setting

E = gp_objfun(logtheta,t,X,covfunc);

% priors
priors = prt_gp_priors(func2str(covfunc));

if iscell(covfunc)
    d = str2double(feval(covfunc{:}));
else
    d = str2double(feval(covfunc));
end
% compute priors
theta = exp(logtheta);
lP  = zeros(d,1);
%dlP = zeros(d,1);
for i = 1:d
    switch priors(i).type
        case 'gauss'
            mu = priors(i).param(1);
            s2 = priors(i).param(2);
            
            lP(i)  = ( -0.5*log(2*pi) - 0.5*log(s2) - 0.5*(theta(i)-mu)^2/s2);
            %dlP(i) = (-(theta(i)-mu) / s2);
            
        case 'gamma'
            a = priors(i).param(1)*priors(i).param(2) + 1;
            b = priors(i).param(2);
            
            lP(i) = (a*log(b) - gammaln(a) + (a - 1)*log(theta(i)) - b*theta(i));
            %%lP(i)  = log(gampdf(theta(i), a, 1/b));
            %dlP(i) =  ((a - 1) / theta(i) - b);
            
        otherwise
            error(['Unknown prior type: ', priors(i).type]);
    end
end

% outputs
nlP = -sum(lP);
%pnlZ  = nlZ + nlP;

E = E + nlP;
end

% -------------------------------------------------------------------------
function [f,F,a] = gp_lap_multiclass(logtheta,covfunc,X,t,f)
% Find mode for Laplace approximation for multi-class classification.
% Derived mostly from Rasmussen & Williams
% Algorithm 3.3 (page 50).
[N,C] = size(t);
if nargin<5, f = zeros(N,C); end;
%if norm(K)>1e8, F=-1e10; return; end

K = covfunc(logtheta,X);

for i=1:32,
    f   = f - repmat(max(f,[],2),1,size(f,2));
    sig = exp(f)+eps;
    sig = sig./repmat(sum(sig,2),1,C);
    E   = zeros(N,N,C);
    for c1=1:C
        D         = sig(:,c1);
        sD        = sqrt(D);
        L         = chol(eye(N) + K.*(sD*sD'));
        E(:,:,c1) = diag(sD)*(L\(L'\diag(sD)));
       %z(c1)     = sum(log(diag(L)));
    end
    M = chol(sum(E,3));

    b = t-sig+sig.*f;
    for c1=1:C,
        for c2=1:C,
            b(:,c1) = b(:,c1) - sig(:,c1).*sig(:,c2).*f(:,c2);
        end
    end

    c   = zeros(size(t));
    for c1=1:C,
        c(:,c1) = E(:,:,c1)*K*b(:,c1);
    end
    tmp = M\(M'\sum(c,2));
    a   = b-c;
    for c1=1:C,
        a(:,c1) = a(:,c1) + E(:,:,c1)*tmp;
    end
    of = f;
    f  = K*a;
   
    %fprintf('%d -> %g %g %g\n', i,-0.5*a(:)'*f(:), t(:)'*f(:), -sum(log(sum(exp(f),2)),1));
    if sum((f(:)-of(:)).^2)<(20*eps)^2*numel(f), break; end
end
if nargout>1
    % Really not sure about sum(z) as being the determinant.
    % hlogdet = sum(z);

    R  = null(ones(1,C));
    sW = sparse([],[],[],N*(C-1),N*(C-1));
    for i=1:N,
        ind         = (0:(C-2))*N+i;
        P           = sig(i,:)';
        D           = diag(P);
        sW(ind,ind) = sqrtm(R'*(D-P*P')*R);
    end
    hlogdet = sum(log(diag(chol(speye(N*(C-1))+sW*kron(eye(C-1),K)*sW))));
    F       = -0.5*a(:)'*f(:) + t(:)'*f(:) - sum(log(sum(exp(f),2)),1) - hlogdet;
    %fprintf('%g %g %g\n', -0.5*a(:)'*f(:) + t(:)'*f(:) - sum(log(sum(exp(f),2)),1), -hlogdet, F);
end
end

% -------------------------------------------------------------------------
function [p Mu SS] = gp_pred_lap_multiclass(logtheta,X,t,covfunc,Xs,Xss,f)
% Predictions for Laplace approximation to multi-class classification.
% Derived mostly from Rasmussen & Williams
% Algorithm 3.4 (page 51).

[N,C] = size(t);

K   = covfunc(logtheta,X);
Ks  = covfunc(logtheta,X,Xs);
kss = covfunc(logtheta,Xss,'diag');

if nargin<7,
    f = gp_lap_multiclass(logtheta,covfunc,X,t);
end

sig = exp(f);
sig = sig./repmat(sum(sig,2)+eps,1,C);
E   = zeros(N,N,C);
for c1=1:C   
    D         = sig(:,c1);
    sD        = sqrt(D);
    L         = chol(eye(N) + K.*(sD*sD'));
    E(:,:,c1) = diag(sD)*(L\(L'\diag(sD) ));
end 
M   = chol(sum(E,3));
try
    os  = RandStream.getGlobalStream;
catch
    os  = RandStream.getDefaultStream;
end
p   = zeros(size(Ks,2),C);
j   = 0;
Mu = zeros(size(Ks,2),C); SS = zeros(C,C,size(Ks,2));
for i=1:size(Ks,2),
    j = j + 1;

    mu = zeros(C,1);
    S  = zeros(C,C);
    for c1=1:C,
        mu(c1) = (t(:,c1)-sig(:,c1))'*Ks(:,i);
        b      = E(:,:,c1)*Ks(:,i);
        c      = (M\(M'\b));
        for c2=1:C,
            S(c1,c2) = Ks(:,i)'*E(:,:,c2)*c;
        end
        S(c1,c1) = S(c1,c1) - b'*Ks(:,i) + kss(i);
    end
    
    % collect latent means and variances
    Mu(i,:)   = mu';
    SS(:,:,i) = S;
    
    s = RandStream.create('mt19937ar','seed',0);
    try
        RandStream.setGlobalStream(s);
    catch
        RandStream.setDefaultStream(s);
    end
    nsamp  = 10000;
    r      = sqrtm(S)*randn(C,nsamp) + repmat(mu,1,nsamp);
    %r      = chol(S)'*randn(C,nsamp) + repmat(mu,1,nsamp);
    % subtract a constant to avoid numerical overflow
    r      = bsxfun(@minus, r, max(r, [], 1));
    r      = exp(r);
    p(j,:) = mean(r./repmat(sum(r,1),C,1),2)';
end
try
    RandStream.setGlobalStream(os);
catch
    RandStream.setDefaultStream(os);
end
end

