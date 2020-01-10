function [y,tw]=nt_mean_over_trials(x,w)
%[y,tw]=nt_mean_over_trials(x,w) - weighted average over trials
%
%  y: weighted average over trials (time*trials)
%  tw: total weight (time*1)
%
%  x: data to average (time*channels*trials)
%  w: weight to apply (time*channels*trials or time*1*trials);

if nargin<2; w=[]; end
if nargin<1; error('!'); end

[m,n,o]=size(x);

if isempty(w); 
    y=mean(x,3); 
    tw=ones(m,n,1)*o;
else
    [mw,nw,ow]=size(w);
    if mw~=m; error('!'); end
    if ow~=o; error('!'); end
    x=nt_unfold(x); 
    w=nt_unfold(w);
    if nw==n;
        x=x.*w;
        x=nt_fold(x,m);
        w=nt_fold(w,m);
        y=sum(x,3)./sum(w,3);
    elseif nw==1;
        x=nt_vecmult(x,w);
        x=nt_fold(x,m);
        w=nt_fold(w,m);
        y=nt_vecmult(sum(x,3),1./sum(w,3));
    end
    tw=sum(w,3);
end

