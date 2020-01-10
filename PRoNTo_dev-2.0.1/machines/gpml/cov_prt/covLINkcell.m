function K = covLINkcell(hyp, A, B, i)

% Linear covariance function with a single hyperparameter. The covariance
% function is parameterized as:
%
% k(x,z) = (x'*z)/t2;
%
% The hyperparameter controls the scaling of the latent function:
%
% hyp = [ log(sqrt(t2)) ]
%__________________________________________________________________________
% This script is based on a script derived from the GPML toolbox, which is
% Copyright (C) Carl Rasmussen and Hannes Nickisch, 2011.

% Written by A Marquand 
% $Id: covLINkernel.m 176 2011-10-20 08:44:21Z amarquan $

if nargin<2, K = '1'; return; end             % report number of parameters
if nargin<3, B = []; end                              % make sure, B exists
xeqz = numel(B)==0; dg = strcmp(B,'diag') && numel(B)>0;   % determine mode

it2 = exp(-2*hyp);                                             % t2 inverse

if iscell(A), A = A{:}; end
if iscell(B), B = B{:}; end

% configure raw kernel from input arguments
if dg % kss
  K = diag(A);
else
  if xeqz % K
    K = A;
  else % Kss
    K = B'; % transpose because gpml expects K(train,test)
  end
end

if nargin<4 % return covariances
  K = it2*K;
else % return derivatives
  if i==1
    K = -2*it2*K;
  else
    error('Unknown hyperparameter')
  end
end
