function [c0,c1,w]=bias_hi_amp(x,exponent)
%[c0,c1]=bias_hi_amp(x,exponent) - covariance with and without hi-amp bias
%
% y: high-amplitude component
% 
% x: data set (2D or 3D)
% exponent: exponent to apply to squared difference (larger=peakier) [default: 0.5]

if nargin<2; exponent=0.5; end

% bias function emphasizes high-amplitude
w=mean((x.^2).^exponent,2);
ww=squeeze(min(w,[],1));
w=nt_vecadd(squeeze(w),-ww');  % subtract minimum from each trial (sharpens the bias function)
w=nt_fold(w(:),size(x,1));

c0=nt_cov(x);
% c1=tscov(vecmult(unfold(x),w));
% save memory:
c1=zeros(size(x,2));
for k=1:size(x,3)
    c1=c1+nt_cov(nt_vecmult(x(:,:,k),w(:,:,k)));
end