function w=nt_find_clipped(x,bounds)
%w=nt_find_clipped_trials(x,bounds) - find clipped trials
%
%  w: 0: clipped, 1: OK
%  
%  x: data (time * channels * trials)
%  bounds: data bounds matrix (nchans x 2) [default: max & min of data in each channel]
% 
% Clipping is presumed to occur on a channel if:
%  - the value is greater than 'bounds' for that channel, or
%  - the value remains at min or max for at least 3 consecutive samples and
%  min and max are not zero
%
% NoiseTools

if nargin<2; bounds=[]; end

[m,n,o]=size(x);

w=zeros(m,1,o);
for k=1:n
    xx=x(:,k,:);
    mx=max(nt_unfold(xx));
    mn=min(nt_unfold(xx));
    ww=zeros(m,1,o);
    ww(find(xx==mx))=1;
    ww(find(xx==mn))=1;
    www=ww(1:end-2,:,:).*ww(2:end-1,:,:).*ww(3:end,:,:);
    if find(www>0) % clipping occured
        w=max(w,ww);
    end
end

w=1-w;