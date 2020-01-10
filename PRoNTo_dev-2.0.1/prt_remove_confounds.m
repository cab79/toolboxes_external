function [Kr, R] = prt_remove_confounds(K,C)

% [Kr, R] = prt_remove_confounds(K,C)
%
% Function to remove confounds from kernel.
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand
% $Id$


n  = size(K,1);
R  = eye(n)-C*pinv(C); %R=eye(n)-C*inv(C'*C)*C'
Kr = R'*K*R;

return