function weights = prt_weights_svm_bin (d,args)
% Run function to compute weights for binary SVM
% FORMAT weights = prt_weights_svm_bin (d,args)
% Inputs:
%       d              - data structure
%           .datamat   - data matrix [Nfeatures x Nexamples]
%           .coeffs    - coefficients vector [Nexamples x 1]
%       args           - function arguments (can be left empty)
% Output:
%       weights        - vector with weights [Nfeatures x 1]
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Mourao-Miranda and M.J.Rosa
% $Id$

SANITYCHECK = true; % turn off for speed

% initial checks
%--------------------------------------------------------------------------
if SANITYCHECK == true
    if isempty(d)
        error('prt_weights_svm_bin:DataEmpty',...
            'Error: ''data'' cannot be empty!');
    else
        if ~isfield(d,'datamat')
            error('prt_weights_svm_bin:noDatamatField',...
                ['Error: ''data'' struct must contain a ''datamat'' '...
                ' field!']);
        end
        if isfield(d,'coeffs')
            if ~isvector(d.coeffs)
                error('prt_weights_svm_bin:CoeffsnoVector',...
                    'Error: ''coeffs'' must be a vector!');
            else
                ncoeffs = length(d.coeffs);
            end
        else
            error('prt_weights_svm_bin:noCoeffsField',...
                ['Error: ''data'' struct must contain ''coeffs'' '...
                ' field!']);
        end
    end
end

% create 1D image
%--------------------------------------------------------------------------
img1d     = zeros(size(d.datamat(1,:)),'single');

% compute weigths
for i=1:ncoeffs   
    
    tmp1 = single(d.datamat(i,:));
    tmp2 = single(d.coeffs(i));
   img1d = img1d + tmp1 * tmp2;
   
end 

% weigths
weights{1} = img1d;