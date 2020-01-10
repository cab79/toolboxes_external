function K_normalised = prt_normalise_kernel(K)

% FORMAT K_normalised = prt_normalise_kernel(K)
%
% This function normalises the kernel matrix such that each entry is 
% divided by the product of the std deviations, i.e.
% K_new(x,y) = K(x,y) / sqrt(var(x)*var(y)) 
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand
% Id: $

d  = diag(K);
K0 = sqrt(repmat(d,[1,size(K,1)]).* repmat(d',[size(K,1),1]));
K_normalised  = K./K0;

end