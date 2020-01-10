function [neglike, neglogprob, prob, mu] = br_regnlike(error, X, y, beta, beta0, s2, tdof)
%BR_REGPROBABILITY computes probabilities of data for linear regression models.
%
%   The input arguments are:
%       error      - error model {'gaussian','laplace','t'}
%       y          - [n x 1] target vector
%       mu         - [n x p] vector of means/locations
%       s2         - [1 x p] residual variance
%       tdof       - [1 x 1] deegrees of freedom (for t-distribution only)
%
%   Return values:
%       neglike    - [1 x p] negative log-likelihood of y for each (beta,beta0,sigma2) sample
%       neglogprob - [n x p] negative log-probability of y | ...
%       prob       - [n x p] probability of y | beta, beta0, s2, tdof
%
%   (c) Copyright Enes Makalic and Daniel F. Schmidt, 2017

mu = bsxfun(@plus, X*beta, beta0);

% If continuous data
if (~strcmp(error,'binomial'))
    e = bsxfun(@minus, mu, y);

    if(strcmp(error,'gaussian'))
        neglogprob = bsxfun(@plus, bsxfun(@rdivide,e.^2,s2*2), (1/2)*log(2*pi*s2));
    elseif(strcmp(error,'laplace'))
        scale = sqrt(s2/2);
        neglogprob = bsxfun(@plus, bsxfun(@rdivide, abs(e), scale), log(2*scale));
    elseif(strcmp(error,'t'))
        nu = tdof;
        neglogprob = -gammaln((nu+1)/2) + gammaln(nu/2) + bsxfun(@plus, (nu+1)/2*log(1 + bsxfun(@rdivide, 1/nu*e.^2, s2)),  log(pi*nu*s2)/2);
    end
    
% Else binary data    
else
    %% numerical constants
    lowerBnd = log(eps); 
    upperBnd = -lowerBnd;
    probLims = [eps, 1-eps];

    %% Compute negative log-likelihood
    mu = constrain(mu, lowerBnd, upperBnd);
    prob = 1./(1 + exp(-mu));

    if any(prob(:) < probLims(1) | probLims(2) < prob(:))
        prob = max(min(mu,probLims(2)),probLims(1));
    end

    neglogprob = -bsxfun(@times, log(prob), y) - bsxfun(@times, log(1.0-prob), (1.0-y));
end    

prob = exp(-neglogprob);
neglike = sum(neglogprob,1);

end

function x = constrain(x,lower,upper)

% Constrain between upper and lower limits, and do not ignore NaN
x(x<lower) = lower;
x(x>upper) = upper;

%% done;
end