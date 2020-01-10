function [images, recons, w, n, mu] = decode_mixtmod_mlr(X, Y, trainidx, testidx, B, Sigma, priordata, T, Pmlr)
%
% Forming reconstructions from the encoding models and prior information
% with Gaussian mixture models and semantic gating
%
% Input
% X: brain data (trials x voxels)
% Y: image data (trials x pixels)
% trainidx: trials used for training
% testidx: trials used for testing
% B: encoding filters
% Sigma: estimated residual variances used to define the model
% priordata: ntrials x pixels x categories prior data
% T: temperature of the mixture
% Pmlr: Probability of belonging to each class in the prior based on other 
% brain area (class x trials)
%
% Output
% images: original images from test set
% recons: reconstruction from brain data of the test set
% w: probability of the category (weight)
% n: MAP reconstruction per category
% mu: means for all components in the prior
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
% Schoenmakers, S., Güçlü, U., Van_gerven, M., & Heskes, T. (2014). Gaussian 
% mixture models and semantic gating improve reconstructions from human 
% brain activity. Name: Frontiers in Computational Neuroscience, 8, 173.
%
% Kindly report any suggestions or corrections to
% s.schoenmakers@donders.ru.nl

if nargin<8,
  T=1;
end

numcat  = size(priordata,3);
nfeatures = size(Y,2); 
nvoxels = size(X,2); 
cat_P = 1/numcat;

ntrain = numel(trainidx);
ntest = numel(testidx);

images = Y(testidx,:);

zet = any(Y);
pixsize = sum(zet);

%remove zeros from prior and Y images
priordata = priordata(:,logical(zet(:)),:);
Y = Y(:,logical(zet(:)));
B = B(logical(zet(:)),:);

% scale everything

[X1,mu1,sigma1] = zscore(X(trainidx,:));
X2 = bsxfun(@rdivide,bsxfun(@minus,X(testidx,:),mu1),sigma1);

[I1,mu2,sigma2] = zscore(Y(trainidx,:));

I2 = bsxfun(@rdivide,bsxfun(@minus,Y(testidx,:),mu2),sigma2);
I2(isnan(I2(:))) = 0; I2(isinf(I2(:))) = 0;

clear X;

mu = zeros(size(priordata,2),numcat);
R = zeros(pixsize,pixsize,numcat);
for i = 1:size(priordata,3)
    
  % calculate the mean per component in the prior
  mu(:,i) = mean(priordata(:,:,i));
  R(1:pixsize,1:pixsize,i) = cov(priordata(:,:,i));

end

% only use voxels with non-zero weights in them
vidx = any(B);
Bt = B(:,vidx);
Xt = X2(:,vidx);

% residual variances for restricted set
Sigma = Sigma(vidx,vidx);

%% reconstruction; we need to compute ni and wi for each trial

BinvS = Bt / Sigma;
D = BinvS * Bt';
zbar = BinvS * Xt';

n = zeros(ntest,pixsize,numcat); % MAP reconstruction per category
w  = zeros(ntest,numcat); % probability of the category (weight)
for i=1:numcat,
    disp(i)

    Uinv = eye(pixsize) + R(:,:,i)*D;
    Umui = Uinv \ mu(:,i);
    Qi = Uinv \ R(:,:,i);
    
    Qizbar = Qi*zbar;
    n(:,:,i) = bsxfun(@plus,Qizbar,Umui)';
    w(:,i) = -sum(log(diag(chol(Uinv*Uinv'))))/2 - mu(:,i)'*D*Umui/2 + ...
        diag(zbar'*Qizbar)/2 + zbar'*Umui +log(Pmlr(i,:)');
    
end

w = bsxfun(@minus,w,max(w,[],2));
    
D2 = cell(1,numel(T)); % vary temperature
for t=1:numel(T)
    if numcat > 1,
      wt = exp(w/T(t));
      wt = bsxfun(@rdivide,wt,sum(wt,2));
      D2{t} = squeeze(sum(bsxfun(@times,shiftdim(n,2),wt'),1));
    else
      D2{t} = n;
    end
end
w = exp(w);
w = bsxfun(@rdivide,w,sum(w,2));

%% rescale everything back to the original

recons = cell(1,numel(T));
for t=1:numel(T)

  recons{t} = bsxfun(@plus,bsxfun(@times,D2{t},sigma2),mu2);

  %% put zero pixels back
  
  tmp = zeros(size(recons{t},1),numel(zet));
  tmp(:,zet) = recons{t};
  recons{t} = tmp;

end

tmp = zeros(size(recons{1},1),numel(zet),numcat);
tmp(:,zet,:) = n;
n = tmp;

tmp = zeros(numel(zet),numcat);
tmp(zet,:) = mu;
mu = tmp;

if numel(recons)==1
  recons = recons{1};
end

end