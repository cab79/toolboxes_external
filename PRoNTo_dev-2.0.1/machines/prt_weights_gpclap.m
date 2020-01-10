function weights = prt_weights_gpclap (d,args)
% Run function to compute weights for linear multiclass classifiers
% FORMAT weights = prt_weights_gpclap (d,args)
% Inputs:
%       d              - data structure
%           .datamat   - data matrix [Nfeatures x Nexamples]
%           .coeffs    - coefficients vector [Nexamples x 1]
%       args           - function arguments (can be empty)
% Output:
%       weights        - vector with weights {Nclass}[Nfeatures x 1]
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M.J.Rosa
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
                [ncoeffs nclass] = size(d.coeffs);
        else
            error('prt_weights_svm_bin:noCoeffsField',...
                ['Error: ''data'' struct must contain ''coeffs'' '...
                ' field!']);
        end
    end
end
weights=cell(nclass,1);
for c = 1:nclass
    % create 1D image
    %----------------------------------------------------------------------
    img1d     = zeros(size(d.datamat(1,:)),'single');
    
    % compute weigths
    for i=1:ncoeffs
        
        tmp1 = single(d.datamat(i,:));
        tmp2 = single(d.coeffs(i,c));
        img1d = img1d + tmp1 * tmp2;
        
    end
    
    % weigths
    weights{c} = img1d;
    
end
