function result = f_johnson_Q(q1,q2,q3,q4)
% - use quantiles to estimate parameters of a Johnson distribution
%
% USAGE: result = f_johnson_Q(q1,q2,q3,q4)
%
% q1 = 1st quantile
% q2 = 2nd quantile
% q3 = 3rd quantile
% q4 = 4th quantile
%
% result = structure of results with the following fields:
%  .coef = parameters as: coef = [gamma delta xi lambda];
%  .type = type of Johnson distribution as: SL, SU, SB, or SN
%
% SEE ALSO: f_johnson_fit, f_johnson_M

% -----Notes:-----
% This function implements Wheeler's (1980's) method for estimating the
% parameters of a Johnson curve using quantiles by calling the 'johnsrnd'
% command in the Matlab Statistics Toolbox. Its main function is to maintain
% consistency with the inputs/outputs of the f_johnson_M and f_johnson_fit
% functions.
% 
% The coefficients of a Johnson curve consist of shape (gamma & delta), 
% location (xi), and scale (lambda) parameters and use the following transforms:
%  SL: Lognormal distribution = exponential transform
%  SU: Unbounded distribution = hyperbolic sine transform
%  SB: Bounded distribution   = logistic transforma
%  SN: Normal distribution    = identity transform
%
% Bounded (SB) curves are bounded on either the lower, upper, or both ends and
% include distributions such as the Beta and Gamma. Unbounded (SU) curves
% include distributions such as the Normal and t distributions.

% -----References:-----
% Johnson, N. L. 1949. Systems of frequency curves generated by methods of
%  translation. Biometrika 36: 149-176.
% Wheeler, R. E. 1980. Quantile estimators of Johnson curve parameters.
%  Biometrika 67(3): 725-728.

% -----Author:-----
% by David L. Jones, Mar-2014
%
% This file is part of the 'JOHNSON CURVE TOOLBOX FOR MATLAB'
% and is released under the BSD 2-clause license.

% -----Set defaults & check input:-----
% Check for scalars:
if (sum([isscalar(q1) isscalar(q2) isscalar(q3) isscalar(q4)])<4)
   error('All inputs must be scalars!');
end
% -------------------------------------

% Estimate parameters based on quantiles:
[~,type,coef] = johnsrnd([q1 q2 q3 q4],0,0);

% Wrap results up into a structure:
result.coef = coef;
result.type = type;
