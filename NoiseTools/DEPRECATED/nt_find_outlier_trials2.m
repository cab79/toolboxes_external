function [idx,d]=nt_find_outlier_trials2(x,criterion,mn,regress_flag)
%[idx,d,mn,idx_unsorted]=nt_find_outlier_trials(x,criterion,mn,regress_flag) - find outlier trials
%
%  idx: indices of trials to keep
%  d: relative deviations from mean
%  
%  x: data (time * channels * trials)
%  criterion: keep trials less than criterion from mean
%  mn: mean (default: calculate from data) 
%  regress_flag: if true regress out mean, rather than subtract
%
%  For example criterion=2 rejects trials that deviate from the mean by
%  more than twice the average deviation from the mean.
%
%  Use nt_find_outlier_trials instead to remove a fixed proportion of trials.
%
%  If no output arguments are specified, plots 'd'.
%

if nargin<2; criterion=inf; end
if nargin<3; mn=[]; end
if nargin<4; regress_flag=0; end
if ndims(x)~=3; error('x should be 3D'); end

if nargout==0;
    [idx,d]=nt_find_outlier_trials2(x,criterion,mn,regress_flag);
    plot(d,'.-');
    xlabel('trial'); ylabel('normalized deviation from mean'); 
    clear idx d mn idx_unsorted
    return
end

[m,n,o]=size(x);
x=reshape(x,m*n,o);

if isempty(mn); mn=mean(x,2); end
if regress_flag
    mn=nt_tsregress(x,mean(x,2));  % regression
else
    mn=repmat(mn(:),1,o);       % mean
end
d=x-mn; % difference from mean
if 0 
    d=sum(d.^2);
else % to save memory
    dd=zeros(1,size(d,2));
    for k=1:size(d,2); dd(k)=sum(d(:,k).^2); end
    d=dd; clear dd;
end
d=d/(sum(x(:).^2)/o);

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


if nargout==0;
    plot(d);
end
