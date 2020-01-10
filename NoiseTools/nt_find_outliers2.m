function w=nt_find_outliers2(x,cutoff,iterations);
%w=nt_find_outliers2(x,cutoff,iterations) - outliers based on mahalanobis distance
%
%  w: mask matrix (0: bad, 1: good)
%
%  x: data 
%  cutoff: outlier if (m.d.)/nchans > cutoff (default: 2)
%  iterations: number of times to iterate (default: 1)
% NoiseTools

if nargin<1; error('!'); return; end
if nargin<2 || isempty(cutoff); cutoff=2; end
if nargin<3 || isempty(iterations); iterations=1; end

[m,n,o]=size(x);
x=nt_unfold(x);

w=ones(size(x,1),1);
if iterations>1
    w=nt_find_outliers2(x(find(w),:),cutoff,iterations-1);
    return
end

d=mahal(x,x(find(w),:));
d=d/n; % normalize by number of channels
w=d<cutoff;

w=nt_fold(w,m);
