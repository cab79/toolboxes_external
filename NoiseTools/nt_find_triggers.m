function tidx=nt_find_triggers(x,threshold,type,guard)
%tidx=nt_find_triggers(x,threshold,type,guard) - find triggers
%
%  tlist: array of trigger indices
%
%  x: trigger channel waveform
%  threshold: trigger threshold
%  type: type of trigger ('up' or 'down',guard);
%  guard: (samples) dead interval to avoid rebound 


if nargin<4 || isempty(guard) ; guard=0; end
if nargin<3 || isempty(type) ; type='up'; end
if nargin<2 threshold=[]; end
if nargin<1; error('!'); end

if strcmp(type,'down');
    tidx=nt_find_triggers(-x,threshold,'up',guard);
    return;
end

if ~strcmp(type,'up'); error('!'); end

if x ~= x(:); error('trigger waveform should be column vector'); end

if isempty(threshold); 
    threshold=(max(x)+min(x))/2;
end

xx=x(1:end-1)<threshold & x(2:end)>=threshold;




tidx=find(xx);
