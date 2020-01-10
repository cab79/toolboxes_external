function nt_widthlines(h,widths)
%nt_widthlines(h,permutation) - apply different widths to lines of plot
%
%  h: handle to plot (default:gca)
%  widths: array of widths to apply to plot
%
% Colors are applied to children of h in reverse order (ie in order of plot
% commands).  May produce unexpected results if there are childern other
% than plot lines.\
% 
% NoiseTools

if nargin<1 || isempty(h); h=gca; end
if nargin<2; error('!'); end

c=get(h,'children');

for k=1:numel(c);
    set(c(numel(c)-k+1),'linewidth', widths(k))
end

