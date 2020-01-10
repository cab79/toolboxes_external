function [tocomps,ii]=nt_peaky(c,x,T,nSmooth)
%[tocomps,ii]=nt_peaky(c,x,T,nSmooth) - find components that maximize peakiness
%
%  tocomps: analysis matrix
%  ii: array of positions
%  
%  c: array of covariance matrices
%  x: data (time*channels or time*channels*trials)
%  T: (samples) time window for each covariance matrix (=dsratio)
%  nSmooth: smoothing factor for covariance matrices
%  
%  This function searches over a series of time intervals for the time interval 
%  that maximizes the ratio between between first and second components of a DSS
%  analysis that contrasts that time interval and the full data.
%  
%  Example:
%    [tocomps]=nt_peaky(c): apply to array of covariance matrices
%    [tocomps]=ny_peaky([],x,T): apply to signal, T is integration window
%    [tocomps]=ny_peaky([],x,T,nSmooth): array is smoothed over nSmoothed samples
nt_greetings;
 
if nargin<1; error('!'); end

topcs=[];
if isempty(c); 
    % we need to calculate the series of covariance matrices
    if nargin<3; error('!'); end
    if nargin<4; nSmooth=1; end
    topcs=nt_pca0(x); % improves condition of matrices
    x=nt_mmat(x,topcs);
    c=nt_xprod(x,'full',T);
    if nSmooth>1;
        c=filter(ones(nSmooth,1),1,c);
        c=c(nSmooth:end,:,:);
    end
end

if ndims(c)~=3 || size(c,2)~=size(c,3)
    error('c has unexpected shape');
end
nComp=size(c,2);

%{
We iterate from 1 to nComp, updating the tocomps matrix at each step.
At each step we iterate over time intervals.  At each time interval we
apply DSS using interval / total power ratio as a bias.  We select the
interval that gives the greatest ratio.

We then recalculate the corresponding DSS solution, put the first component
into the tocomps matrix, and recurse on the remainder.

Calculations are done on covariance matrices rather than actual data.
%}

c0=squeeze(mean(c));
if ~isempty(topcs)
    tocomps=topcs;
else
    tocomps=eye(nComp);
end
iComp=1;
[m n o]= size(c); % c is ntime x nComp x nComp
while iComp<nComp
    
    % find the index that gives the smallest D2/D1 ratio
    ratio=zeros(1,size(c,1));
    for iIter=1: size(c,1);
        if size(c0,1)>150   % test to optimize speed
            N=10;
            [V,D]=eigs(squeeze(c(iIter,:,:)),c0,N);
        else
            [V,D]=eig(squeeze(c(iIter,:,:)),c0);
        end
        D=sort(diag(D),'descend');
        if 0
            ratio(iIter)=D(2)/D(1);
        else
            ratio(iIter)=mean(D(2:end))/D(1);
        end
    end
    [~,idx]=min(ratio);
    ii(iComp)=idx;
    
    %plot(ratio); pause
    
    %[cond(c0) cond(squeeze(c(idx,:,:)))]
    
    % DSS
    [todss,pwr0,pwr1]=nt_dss0(c0,squeeze(c(idx,:,:)),[],10.^-10);
    profile=pwr1./pwr0;
    score(iComp)=profile(1)/profile(2);
    
    % update tocomps
    tocomps=[tocomps(:,1:iComp-1), tocomps(:,iComp:end)*todss];
    nComp=size(tocomps,2);
    
    % update covariance matrices
    if iComp<nComp
        todss=todss(:,2:end);
        c0=todss'*c0*todss;
        %{
        For speed, we reshape c to apply left and right multiplication to all
        matrices at once.
        %}
        [m n o]=size(c);
        c=reshape(permute(c,[2 3 1]), [n o*m]); % concatenate horizontally
        c=todss'*c; % left multiply
        n=size(c,1);
        c=reshape(permute(reshape(c,[n o m]),[1 3 2]),[m*n o]); % concatenate vertically
        c=c*todss; % right multiply
        o=size(c,2);
        c=permute(reshape(c, [n m o]),[2 1 3]) ; % reshape to ntime x nComp x nComp
    end
    iComp=iComp+1;
end

        
        
        