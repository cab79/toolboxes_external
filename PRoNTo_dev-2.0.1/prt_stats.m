function stats = prt_stats(model, tte, nk)
% Function to compute predictions machine performance statistcs statistics
%
% Inputs:
% ----------------
% model.predictions: predictions derived from the predictive model
% model.type:        what type of prediction machine (e.g. 'classifier','regression')
%
% tte: true targets (test set)
% nk:  number of classes if classification (empty otherwise)
% flag:  'fold' for statistics in each fold
%         'model' for statistics in each model
% 
% Outputs:
%-------------------
% Classification:
% stats.con_mat: Confusion matrix (nClasses x nClasses matrix, pred x true)
% stats.acc:     Accuracy (scalar)
% stats.b_acc:   Balanced accuracy (nClasses x 1 vector)
% stats.c_acc:   Accuracy by class (nClasses x 1 vector)
% stats.c_pv:    Predictive value for each class (nClasses x 1 vector)
%
% Regression:
% stats.mse:     Mean square error between test and prediction
% stats.corr:    Correlation between test and prediction
% stats.r2:      Squared correlation
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand
% $Id$

% FIXME: is any code using the 'flags' input argument? 
if ~isfield(model,'type')
    warning('prt_stats:modelDoesNotProvideTypeField',...
        'model.type not specified, defaulting to classifier');
    model.type = 'classifier';
end

switch model.type
    case 'classifier'
        
        stats = compute_stats_classifier(model, tte, nk);
        
    case 'regression'
        
        stats = compute_stats_regression(model, tte);
        
    otherwise
        error('prt_stats:unknownTypeSpecified',...
            ['No method exists for processing machine: ',machine.type]);
end

end

% -------------------------------------------------------------------------
% Private functions
% -------------------------------------------------------------------------

function stats = compute_stats_classifier(model, tte, k)

k = max(unique(k));        % number of classes

stats.con_mat = zeros(k,k);
for i = 1:length(tte)
    true_lb = tte(i);
    pred_lb = model.predictions(i);
    stats.con_mat(pred_lb,true_lb) = stats.con_mat(pred_lb,true_lb) + 1;
end

Cc = diag(stats.con_mat);   % correct predictions for each class
Zc = sum(stats.con_mat)';   % total samples for each class (cols)
nz = Zc ~= 0;               % classes with nonzero totals (cols)
Zcr = sum(stats.con_mat,2); % total predictions for each class (rows)
nzr = Zcr ~= 0;             % classes with nonzero totals (rows)

stats.acc       = sum(Cc) ./ sum(Zc);
stats.c_acc     = zeros(k,1);
stats.c_acc(nz) = Cc(nz) ./ Zc(nz);
stats.b_acc     = mean(stats.c_acc);
stats.c_pv      = zeros(k,1);
stats.c_pv(nzr) = Cc(nzr) ./ Zcr(nzr); 

% confidence interval
% TODO: check IID assumption here (chunks in run_permutation.m)
% before applying tests, and give nans if not applicable...
[lb,ub] = computeWilsonBinomialCI(sum(Cc),sum(Zc));
stats.acc_lb=lb;
stats.acc_ub=ub;
end

function stats = compute_stats_regression(model, tte)

if numel(tte)<3
    stats.corr = NaN;
    stats.r2 = NaN;
else
    coef = corrcoef(model.predictions,tte);
    stats.corr = coef(1,2);
    stats.r2 = coef(1,2).^2;
end
stats.mse  = mean((model.predictions-tte).^2);
stats.nmse = mean((model.predictions-tte).^2)/(max(tte)-min(tte));
end

function [lb,ub] = computeWilsonBinomialCI(k,n)
% Compute upper and lower 5% confidence interval bounds
% for a binomial distribution using Wilson's 'score interval'
%
% IN
%   k: scalar, number of successes
%   n: scalar, number of samples
%
% OUT
%   lb: lower bound of confidence interval
%   ub: upper bound of confidence interval
%
% REFERENCES
% Brown, Lawrence D., Cai, T. Tony, Dasgupta, Anirban, 1999.
%  Interval estimation for a binomial proportion. Stat. Sci. 16, 101?133.
% Edwin B. Wilson, Probable Inference, the Law of Succession, and
%   Statistical Inference, Journal of the American Statistical Association,
%   Vol. 22, No. 158 (Jun., 1927), pp. 209-212

alpha=0.05;

l=spm_invNcdf(1-alpha/2,0,1); % 
p=k/n;                    % sample proportion of success
q=1-p;

% compute terms of formula
firstTerm=(k+(l^2)/2)/(n+l^2);
secondTerm=((l*sqrt(n))/(n+l^2))*sqrt(p*q+((l^2)/(4*n)));

% compute upper and lower bounds
lb=firstTerm-secondTerm;
ub=firstTerm+secondTerm;

end
