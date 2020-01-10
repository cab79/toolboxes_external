function [c0,c1]=bias_mask(x,mask)
%[c0,c1]=bias_hi_amp(x,mask) - covariance of masked signal
%
% c0: covariance of part for which mask<0
% c1: covariance of part for which mask>0
% 
% x: data set (time X channels or time X channels X trials)
% mask: mask function (time X 1 or time X 1 X trials)

if nargin<2; error('!'); end

x=nt_unfold(x);
mask=nt_unfold(mask);


c0=nt_cov(nt_vecmult(x(find(mask<=0),:),mask(find(mask<=0))));
c1=nt_cov(nt_vecmult(x(find(mask>0),:),mask(find(mask>0))));