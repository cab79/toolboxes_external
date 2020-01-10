function w = prt_KRR(K,t,reg)

% w = prt_KRR(K,t,reg)
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Carlton Chu
% $Id$

w = (K+reg*eye(size(K)))\t;
return