function [idx,d,cc]=nt_find_outlier_trials3(x,criterion,norm_flag)
%[idx,d]=nt_find_outlier_trials3(x,criterion,norm_flag) - find outlier trials
%
%  idx: indices of trials to keep
%  d: relative deviations from mean
%  
%  x: data (time * channels * trials)
%  criterion: keep trials less than criterion from mean
%  norm_flag: if true divide each frame by its RMS [default:0]
%
%  This version compares covariance matrices of sample-normalized data.
%
%  If no output arguments are specified, plots 'd'.
%

if nargin<3||isempty(norm_flag); norm_flag=0; end
if nargin<2||isempty(criterion); criterion=inf; end
if ndims(x)~=3; error('x should be 3D'); end

if nargout==0;
    [idx,d]=nt_find_outlier_trials3(x,criterion,norm_flag);
    plot(d, '.-');
    xlabel('trial'); ylabel('normalized deviation from mean'); 
    clear idx d mn idx_unsorted
    return
end

[m,n,o]=size(x);
x=nt_unfold(x);
if norm_flag
    x=bsxfun(@times,x,1./(eps+sqrt(mean(x.^2,2)))); % normalize by dividing each sample by rms
end
x=nt_fold(x,m);

cc=zeros(n,n,o);
for iTrial=1:o
    cc(:,:,iTrial)=nt_cov(x(:,:,iTrial));
end

cc=reshape(cc,n*n,o);
cc=bsxfun(@times,cc,1./(eps+sqrt(mean(cc.^2,1)))); % normalize by dividing each trial cv by rms

d=sqrt(mean( bsxfun(@minus,cc,mean(cc,2) ).^2 )) ./  sqrt(mean(mean(cc,2).^2));

% 
% 
% if isempty(mn); mn=mean(x,2); end
% if isnan(mn)
%     mn=nt_tsregress(x,mean(x,2));  % distance from regression
% else
%     mn2=repmat(mn(:),1,o);       % distance from mean
% end
% d=x-mn2;


%d=mean(d.^2);
%d=d/sqrt(mean(d.^2));
%d=mean(d.^2)/mean(mn.^2);

idx=find(d<criterion);

% [dd,idx]=sort(d,'ascend');
% idx=idx(find(dd<criterion));
% idx_unsorted=idx;
% idx=sort(idx); % put them back in natural order
% mn=mean(x(:,idx),2);


