function [idx,n_clipped]=nt_find_clipped_trials(x,threshold,bounds)

error('not yet implemented')

%[idx,n_clipped]=nt_find_clipped_trials(x,threshold,bounds) - find clipped trials
%
%  idx: indices of trials to keep
%  n_clipped: number of clipped samples per trial
%  
%  x: data (time * channels * trials)
%  threshold: admissible proportion of clipped samples [default: 0]
%  bounds: clipping bounds [default: max & min of data]
% 
%
% NoiseTools

if nargin<2; threshold=0; end
if nargin<3; bounds=[]; end
if threshold>1; error('proportion should be within [0 1])'); end

if nargout==0;
    [idx,n]=nt_find_outlier_trials(x,proportion,mn);
    plot(n);
    xlabel('trial'); ylabel('number clipped'); 
    clear idx n
    return
end

[m,n,o]=size(x);

mx=max(nt_unfold(x));
mn=min(nt_unfold(x));

x1=repmat(mx,[m,1,o])-x; % zero for clipped, pos for OK
x1=min(x1,[],2); % min over channels
x2=x-repmat(mn,[m,1,o]); % zero for clipped, pos for OK
x2=min(x2,[],2); % min over channels
xx=min(x1,x2); % zero for clipped samples

y=zeros([m,1,o]);
y(find(xx==0))=1;
y=reshape(y,[m,o]); y=sum(y);

plot(y);

