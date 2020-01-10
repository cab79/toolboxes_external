function [orig, recons] = decode_unimod(X, Y, trainidx, testidx, B, Sigma, R)
%
% Forming reconstructions from the encoding model and prior information
%
% Input
% X: brain data (trials x voxels)
% Y: image data (trials x predictors)
% trainidx: trials used for training
% testidx: trials used for testing
% B: encoding filters (predictors x voxels)
% Sigma: estimated residual variances used to define the model (voxels x voxels)
% R: prior covariance matrix (predictors x predictors)
%
% Output
% orig: original images from test set
% recons: reconstruction from brain data of the test set
%
% ------------------------------------------------------------------------
% regularized linear regression toolbox, December 2014
% (c) Sanne Schoenmakers, Tom Heskes, Marcel van Gerven
% Donders Institute for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is hereby
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%
% Please refer to the following paper:
% Schoenmakers, S., Barth, M., Heskes, T., & van Gerven, M. (2013). Linear 
% reconstruction of perceived images from human brain activity. NeuroImage,
% 83, 951-961.
%
% Kindly report any suggestions or corrections to
% s.schoenmakers@donders.ru.nl


nfeatures = size(Y,2);
nvoxels = size(X,2);

% scale everything
% [X1,mu1,sigma1] = zscore(X(trainidx,:));
% X2 = bsxfun(@rdivide,bsxfun(@minus,X(testidx,:),mu1),sigma1); % correct test data to mean and std of training data
% [I1,mu2,sigma2] = zscore(Y(trainidx,:));
% I2 = bsxfun(@rdivide,bsxfun(@minus,Y(testidx,:),mu2),sigma2);
X2 = X(testidx,:); % correct test data to mean and std of training data
I2 = Y(testidx,:);

clear X Y

% DECODING

% only use voxels with non-zero weights in them
vidx = any(B,1);
Bt = B(:,vidx);

% residual variances for restricted set
Sigma = Sigma(vidx,vidx);

% Compute covariance 
BS = Bt/Sigma;
if nfeatures < nvoxels
  
  Q = inv(inv(R) + BS*Bt');
  % BS*Bt' is predictors x predictors (same as R)
  
else
  
  % using matrix inversion lemma
  Q = R - R*Bt*inv(Sigma + Bt'*R*Bt)*Bt'*R;

end

M1 = Q*BS; 

D2 = zeros(size(X2,1),size(Bt,1));
for j=1:size(X2,1)
  D2(j,:) = M1*X2(j,vidx)';
end
 
%% rescale everything back to the original

% orig = bsxfun(@plus,bsxfun(@times,I2,sigma2),mu2);
% recons = bsxfun(@plus,bsxfun(@times,D2,sigma2),mu2);
orig = I2;
recons = D2;

end