function y=nt_hi_amp(x,exponent)
%y=nt_hi_amp(x) - decompose into high-amplitude components


if nargin<1; error('!'); end
if nargin<2 || isempty(exponent); exponent=0.5; end

if ndims(x)==3; 
    y=nt_fold(nt_hi_amp(nt_unfold(x)),size(x,1)); 
    return
end

y=zeros(size(x));
for k=1:size(x,2)-1
    x=nt_normcol(x);
    [c0,c1]=nt_bias_hi_amp(x,exponent);
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    
    figure(100); clf; 
    plot(pwr1./pwr0,'.-'); 
    ylim([0 inf]); 
    xlabel('component'); ylabel('score'); title ('DSS for high samplitude');
    
    x=nt_mmat(x,todss);
    y(:,k)=x(:,1);
    x=x(:,2:end);
    if size(x,2)==0; 
        break;
    end
end

y=y(:,1:k+1);