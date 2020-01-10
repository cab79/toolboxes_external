function [mn,sd,all]=nt_bsmean(x,N)
%[mn,sd,all]=nt_bsmean(x,N) - calculate mean, estimate sd using bootstrap
%
%  mn: mean of x over second dimension
%  sd: standard deviation from mn of bootstrap trials
%  all: matrix of all bootstrap trials
%  
%  x: matrix of observations (time X repetitions)
%  N: number of bootstrap trials [default: 100]

if nargin <2; N=100; end

if ndims(x)>2; 
    x=squeeze(x);
    if ndims(x)>2; 
        error('data must be at most 2D'); 
    end
end

[m,n]=size(x);
all=zeros(m,N);
for k=1:N
    idx=ceil(n*rand(1,n));
    all(:,k)=mean(x(:,idx),2);
end

mn=mean(x,2);
sd=sqrt(mean((all-repmat(mn,1,N)).^2,2));


