function [todss,pwr0,pwr1]=nt_dss1(x,w,keep1,keep2,sns_flag)
%[todss,pwr0,pwr1]=nt_dss1(x,w,keep1,keep2,sns_flag) - evoked-biased DSS denoising
%
%  todss: denoising matrix
%  pwr0: power per component (raw)
%  pwr1: power per component (averaged)
%
%  x: data to denoise (time * channels * trials)
%  w: weight
%  keep1: (in DSS0) number of PCs to retain (default: all)
%  keep2: (in DSS0) ignore PCs smaller than keep2 (default: 10.^-12)
%  sns_flag: if true, apply sns to average
%
%  The data mean is NOT removed prior to processing.
%
% NoiseTools

if nargin<5; sns_flag=[]; end
if nargin<4; keep2=10.^-12; end
if nargin<3; keep1=[]; end
if nargin<2; w=[]; end
if nargin<1; error('!'); end

if ndims(x)<3; error('x should be 3D'); end

[m,n,o]=size(x);
%[x,mn]=nt_demean(x,w);                            % remove weighted mean    

if isempty(w)% weighted average over trials (--> bias function for DSS)
    [c0,nc0]=nt_cov(x);
    c0=c0/nc0;
    [c1,nc1]=nt_cov(mean(x,3)); 
    c1=c1/nc1;
else
    % weighted average over trials (--> bias function for DSS)
    if 1
        [xx,ww]=nt_mean_over_trials(x,w);
        if ~isempty(sns_flag); xx=nt_sns(xx,10,[],w); end
        ww=min(ww,[],2);
    else
        xx=mean(x,3);
        ww=ones([m,1]);
    end
    % covariance of raw and biased data
    [c0,nc0]=nt_cov(x,[],w);
    c0=c0/nc0;
    [c1,nc1]=nt_cov(xx,[],ww); 
    c1=c1/nc1;
end

% derive DSS matrix
[todss,pwr0,pwr1]=nt_dss0(c0,c1,keep1,keep2);



