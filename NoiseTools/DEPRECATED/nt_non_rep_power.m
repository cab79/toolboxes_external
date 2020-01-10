function [r,pwr0,pwr1]=nt_non_rep_power(x,nIterations)
%[r,pwr0,pwr1]=nt_non_rep_power(x,nIterations,nComps) - find non repeating components
%
%  r: component matrix (time*comp*trial)
%  pwr0, pwr1: from DSS
%
%  x: data (time*comp*trial)
%  nIterations: number of times to iterate DSS [default: 2]
%  nComps: number of components to return [default: all]
%

if nargin<3||isempty(nComps); nComps=inf; end
if nargin<2||isempty(nIterations); nIterations=2; end

if ~iscell(x) ; x={x}; end



for k=1:nIterations
    [m,n,o]=size(x{1});
    c0=zeros(n);
    c1=zeros(n);
    for k=1:numel(x)
        x{k}=x{k}/sqrt(mean(x{k}(:).^2)); % divide by total power over this condition
        c0=c0+nt_cov(nt_normpage(x{k}));
        c1=c1+nt_cov(x{k});
    end
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    for k=1:numel(x)
        x{k}=nt_mmat(x{k},todss);
    end
end

r=x;   