function score=nt_repeatability(x,demean_flag)
% [score]=nt_repeatability(x,demean_flag) - repeatability score
%
%  score: ratio of power of mean to total power
%
%  x: data (time * channels * repeats)
%  demean_flag: 0: no demean, 1: demean, 2: remove mean
%  of e
%
%  If no out argument is given, plot result.

if nargin<1; error('!'); end
if ndims(x)~=3; error('data must be 3D'); end

if nargin<2 || isempty(demean_flag)
    demean_flag=1;
end

if demean_flag==1
    x=nt_demean(x);
elseif demean_flag==2
    x=nt_demean2(x);
end

score=mean(mean(x,3).^2)./mean(nt_unfold(x.^2));

if nargout == 0 ;
    plot(score); xlabel('component'); ylabel('score'); title('repeatability');
    score=[];
end

