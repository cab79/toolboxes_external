function [rms,sd,all]=nt_bsrms(x,N)
%[rms,sd,all]=nt_bsrms(x,N) - calculate rms, estimate sd using bootstrap
%
%  rms: rms over second dimension of x
%  sd: standard deviation of rms calculated by bootstrap
%  all: matrix of all trials
%  
%  x: matrix of observations (time X repetitions)
%  N: number of bootstrap trials [default: 100]

if nargin <2; N=100; end

if ndims(x)>2; error('data must be at most 2D'); end

[m,n]=size(x);
all=zeros(m,N);
for k=1:N
    idx=ceil(n*rand(1,n));
    all(:,k)=sqrt(mean(x(:,idx).^2,2));
end

rms=sqrt(mean(x.^2,2));
sd=sqrt(mean((all-repmat(rms,1,N)).^2,2));


