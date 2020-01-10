function nt_bsplot(x,band,style,abscissa,zeroflag,rmsflag)
%nt_bsplot(x,sds,style,abscissa,zeroflag,rmsflag) - plot average with bootstrap standard deviation
%
%  x: data to plot (time * trials, or time * 1 * trials)
%  band: width of band to plot in standard deviations (default: 2)
%  style: 'zerobased' (default) or 'meanbased'
%  abscissa: use this vector as plot's abscissa (as in 'plot(abscissa,x)' )
%  zeroflag: if 1 draw zero line (default: 1)
%  rmsflag: if 1 use RMS instead of mean (default==0)
%
%  Bootstrap uses N=1000 iterations.
% 
%  Example:
%    nt_bsplot(x)
%  where x is time*trials will plot the average of x over trials, together
%  with +/- 2SDs of the bootstrap resampling.
%
% NoiseTools.

if nargin<6 || isempty(rmsflag) ; rmsflag=0; end
if nargin<5 || isempty(zeroflag) ; zeroflag=1; end
if nargin<4; abscissa=[]; end
if nargin<3 || isempty(style); style='zerobased'; end
if nargin<2 || isempty(band); band=2; end

x=squeeze(x);
if ndims(x)>2; error('X should have at most 2 non-singleton dimensions'); end
[m,n]=size(x);
if n<2; error('bootstrap resampling requires more than 1 column'); end
if isempty(abscissa); abscissa=1:m; end
if numel(abscissa) ~= size(x,1); error('abscissa should be same size as x'); end

if rmsflag
    [a,b]=nt_bsrms(x);
else
    N=1000;
    [a,b]=nt_bsmean(x,N);
end
b=b*band;


if strcmp(style,'zerobased');
    Y=[b;-flipud(b)]';
elseif strcmp(style,'meanbased');
    Y=[b+a;flipud(-b+a)]';
else
    error('!');
end
abscissa=abscissa(:);
X=[abscissa;flipud(abscissa)];
C=0.7*[1 1 1];
fill(X,Y,C,'edgecolor','none');
hold on;
plot(abscissa,a*0,'k');
plot(abscissa,a, 'b'); 
hold off

% return
% 
% abscissa2=linspace(min(abscissa),max(abscissa),m*2);
% plot(abscissa2,b,'g');
% c=get(gca,'children'); set(c(1),'color',[.7 .7 .7])
% hold on;
% plot(abscissa,a,'b');
% if zeroflag; 
%     plot(abscissa,0*a,'k'); 
%     c=get(gca,'children'); set(c(1),'color',[.5 .5 .5])
% end
% %set(gca,'xlim',[1 m])
% hold off
