function [todss,K]=nt_kurtosis(x,nIterations,exponent,w,smooth)
%[todss,K]=nt_kurtosis(x,nIterations,exponent,w,smooth)- find high-kurtosis components
%
% x: data
% nIterations: number of times to iterate calculation (default: 5)
% exponent: exponent to apply to mask [default: 1]
% w: weighting function
% smooth: samples, smoothing to apply (twice) to mask

if nargin<1; error('!'); end
if nargin<2||isempty(nIterations); nIterations=5; end
if nargin<3||isempty(exponent); exponent=2; end
if nargin<4; w=[]; end  
if nargin<5; smooth=[];end

if ndims(x)>2; x=nt_unfold(x); end
 
KK={};
K0=kurtosis(x);
if nIterations==1
    
    %{
    We choose the channel with highest kurtosis to define a mask to
    emphasize high amplitude intervals. DSS then finds components with
    variance maximal within those intervals.
    %}
    
    K=kurtosis(x);
    [~,idx]=max(K);
    mask=x(:,idx).^exponent;
    mask=filtfilt(ones(5,1),1,mask);
    if ~isempty(smooth)&&smooth>1; 
        mask=filtfilt(ones(1,smooth),1,mask); % filtfilt keeps mask aligned with data
    end
    c0=nt_cov(x,[],w);
    c1=nt_cov(nt_demean(x.*bsxfun(@times,x,mask),w),[],w);
    [todss]=nt_dss0(c0,c1);
    KK{1}=kurtosis(nt_mmat(x,todss));
    
else
    
    %{
    Apply DSS repeatedly, concatenating the DSS matrices. We hope (no
    guarantee) that the solution will reveal components with high kurtosis.
    %}
    
    z=x;
    todss=eye(size(x,2));
    for k=1:nIterations
        [T]=nt_kurtosis(z,1,exponent,w,smooth);
        z=nt_mmat(z,T);
        todss=todss*T;
        KK{k}=kurtosis(z);
    end
    
end
K=KK{end};

if nargout==0; 
    % don't return anything, just plot
    figure(101); clf
    subplot 121;
    plot(K0); hold on
    for k=1:nIterations
        plot(KK{k}); 
    end
    nt_colorlines;
    title('kurtosis score for each iteration'); xlabel('component');
    legend(num2str((0:nIterations)')); set(gca,'yscale','log')
    subplot 222;
    z=nt_mmat(x,todss(:,1));
    plot(z); title('first DSS component'); 
    subplot 224;
    z=nt_mmat(x,todss(:,end));
    plot(z); title('last DSS component'); 
    clear c0 c1 todss pwr0 pwr1 K
end

