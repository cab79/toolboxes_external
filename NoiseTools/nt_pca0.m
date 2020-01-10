function [topcs,pwr]=nt_pca0(x,shifts,nkeep,threshold,w)
%[topcs,pwr]=nt_pca0(x,shifts,nkeep,threshold,w) - time-shift pca
%
%  topcs: matrix to convert data to PCs
%  pwr: power per PC
%
%  x: data matrix
%  shifts: array of shifts to apply
%  nkeep: number of PCs to keep
%  w: weight (see nt_cov)
%  threshold: remove components with normalized eigenvalues smaller than threshold (default: 0)
%
% mean is NOT removed prior to processing


if nargin<1; error('!'); end
if nargin<2||isempty(shifts); shifts=[0]; end
if nargin<3; nkeep=[]; end
if nargin<4||isempty(threshold); threshold=0; end
if nargin<5; w=[]; end

[m,n,o]=size(x);

% remove mean
%x=fold(demean(unfold(x)),size(x,1));

% covariance
if isempty(w);
    c=nt_cov(x,shifts);
else
    c=nt_cov(x,shifts,w);
end

% PCA matrix
if ~isempty(nkeep)
    topcs=nt_pcarot(c,nkeep);
else
    topcs=nt_pcarot(c);
end

%if ~isempty(nkeep); topcs=topcs(:,1:nkeep); end

% power per PC
pwr=diag(topcs'*c*topcs)/(m*o);
idx=find(pwr/max(pwr)>threshold);
pwr=pwr(idx);
topcs=topcs(:,idx);

% matrix to normalized PCs
%topcs=topcs*diag(1./sqrt(pwr));
