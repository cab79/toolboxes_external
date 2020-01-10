function weights = prt_weights_sMKL_reg(d,args)
% Run function to compute weights for binary MKL
% FORMAT weights = prt_weights_sMKL (d,args)
% Inputs:
%       d               - data structure
%           .datamat    - data matrix [Nfeatures x Nexamples]
%           .coeffs     - coefficients vector [Nexamples x 1]
%           .betas      - kernel weights
%           .idfeat_img - cell with indece
%       args            - function arguments (can be left empty)
% Output:
%       weights         - vector with weights [Nfeatures x 1]
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Mourao-Miranda

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


% compute weigths

img1d     = zeros(size(d.datamat(1,:)),'single');

for k=1:length(args.betas)
    
    betas = single(args.betas(k));
    index_k = args.idfeat_img{k};
    
    if ~isempty(index_k) && betas~=0
        if ~isfield(args,'flag') || ~args.flag
            for i=1:ncoeffs
            
            tmp1 = single(d.datamat(i,index_k));
            tmp2 = single(d.coeffs(i));
            
            img1d(index_k) = img1d(index_k) + tmp1 * tmp2;
            
            end   
        else
            img1d(index_k) = ones(1,length(index_k));
        end        
        
        img1d(index_k) = betas * img1d(index_k);
    end
end

% weigths
weights{1}  = img1d; % originally, it was: weights  = img1d