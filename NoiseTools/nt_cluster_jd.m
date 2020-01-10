function [IDX,TODSS,SCORE]=nt_cluster_jd(x,dsr,smooth,flags,init,verbose)
%[IDX,todss,SCORE]=nt_cluster_jd(x,dsr,smooth,flags,init,verbose) - cluster with joint diagonalization
%
%  IDX: cluster ownership (IDX{1}: low amp, IDX{2{: high amp)
%  TODSS: DSS matrix (1st column --> discriminating component)
%  SCORE: score (smaller means better contrast)
%
%  x: data (time*channels)
%  dsr: downsample ratio for cross product series
%  smooth: further smooth cross-product series
%  flags: 'norm', 'norm2': give each dsr-sized slice the same weight
%  init: provide initial clustering
%  verbose: display & plot (default=no)
%
% See nt_bias_cluster, nt_cluster1D


if nargin<2; error('!'); end
if nargin<3 ||isempty(smooth); smooth=1; end
if nargin<4 ||isempty(flags); flags=[]; end
if nargin<5; init=[]; end
if nargin<6||isempty(verbose); verbose=0; end

if ~nargout; disp('entering nt_cluster_jd...'); end

if ndims(x)>2 || size(x,2) ==1;
    error('should be 2D matrix');
end


%{
 Calculate the time series of cross products (terms of the covariance matrix).
 This time series has coarser temporal resolution than x by a factor dsr.
%}
[xx,ind]=nt_xprod(x,'lower',dsr);  

% figure(2); clf;
% subplot 211;
% plot(xx)

% option: give each slice the same weight (counters amplitude variations)
if find(strcmp(flags,'norm'))
    xx=nt_normrow(xx);
end
if find(strcmp(flags,'norm2'))
    xx=norm2(xx,size(x,2),ind);
end

% subplot 212; 
% plot(xx); 
% pause;

smooth
size(xx)
xx=nt_smooth(xx,smooth,[],1);

%{
Cluster each column the time series of cross products, 
choose the column with best score (reduction in energy), 
and use it's cluster index to initialize the first JD analysis.
%}
if isempty(init)
    [C,A,score]=nt_cluster1D(xx);
    [~,idx]=min(score); % select column with best score (tightest clusters)
    A=A(:,idx); 
        
    % upsample the cluster ownership index so we can apply it to x
    A=repmat(A',[dsr,1]);
    A=A(:);
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==0);
else
    IDX{1}=init;
end

% initial DSS to contrast clusters
c0=nt_cov(x);
c1=nt_cov(x(IDX{1},:));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
todss2=todss(:,[1 end]); % keep only first and last components
z=nt_mmat(x,todss2);

PLOT_FIG2=0;
if PLOT_FIG2
    figure(2);  clf; set(gcf, 'name','in nt_cluster_jd');
    A=zeros(size(x,1),1); A(IDX{1})=1;
    subplot 411; plot(x); title('data');
    subplot 412; plot(A,'.-'); title('initial cluster map');
    subplot 413; plot(z(:,1)); title('initial DSS1');
    subplot 414; plot(z(:,2)); title('initial DSS2');
    drawnow; pause;
end

% iterate until stable
old_IDX=IDX{1};
for k=1:10

    [zz,ind]=nt_xprod(z,[],dsr);
    zz=zz(:,1:2);       % keep only the squares
    zz=log2(zz+eps);    % log to make it sensitive to power ratio
    [C,A]=nt_cluster1D(zz);
    [~,idx]= max(diff(C));    % the best component is the one with greatest ratio
    A=A(:,idx);
    
    %plot(nt_smooth(A,smooth, [],1)); 
    
    A=double(nt_smooth(A,smooth, [],1)>=1/smooth); % extend ownership to include effect of smoothing

    % upsample the cluster ownership index so we can apply it to x
    A=repmat(A',[dsr,1]); % upsample 
    A=A(:); 
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==0); % 0: low values, 1: high values
    
    
    % DSS to contrast clusters
    c0=nt_cov(x)/size(x,1);
    c1=nt_cov(x(IDX{1},:))/size(x(IDX{1},:),1);
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    z=nt_mmat(x,todss(:,[1 end])); % keep first and last

    if ~nargout||verbose; 
        disp(['low amp cluster: ', num2str((100*numel(IDX{1})/size(x,1)), 2), ' % of samples, power ratio: ' num2str(pwr1(end)/pwr0(end), 3)]); 
    end

    if PLOT_FIG2
        figure(2);  
        subplot 515; plot(A,'.-'); title('final cluster map');
    end
    if all(size(old_IDX)==size(IDX{1})) && all(old_IDX==IDX{1}); break; end
    old_IDX=IDX{1};
end 
IDX{2}=setdiff(1:size(x,1), IDX{1});

% score
[TODSS,pwr0,pwr1]=nt_dss0(c0,c1);
SCORE=pwr1(end)/pwr0(end);

nargout

if nargout==0||verbose;
    
    % no output, just plot

    z1=nt_mmat(x,TODSS(:,1));
    z2=nt_mmat(x,TODSS(:,end));

    figure(101); clf ;
    subplot 221;
    plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title('DSS cluster [low amp] vs all');
    subplot 222;
    wsize=min(1024,size(z2,1));
    %nt_spect_plot(z1/sqrt(mean(z1(:).^2)),wsize,[],[],1);
    hold on
    nt_spect_plot(z2/sqrt(mean(z2(:).^2)),wsize,[],[],1);
    nt_spect_plot(x/sqrt(mean(x(:).^2)),wsize,[],[],1);
    xlim([0 .5]);
    nt_linecolors([],[1 3 2]);
    legend('last','all'); legend boxoff
    hold off

    z=nt_mmat(x,todss); 
    z=nt_normcol(z);
    subplot 223; imagescc(nt_cov(z(IDX{1},:))); title('cluster [low amp]'); 
    subplot 224; imagescc(nt_cov(z)-nt_cov(z(IDX{1},:))); title('cluster [high amp]');

    
    figure(102); clf
    if 0
        subplot 311;
        plot(x); hold on
        xx=x; xx(IDX{1},:)=nan;
        plot(xx,'k');
        axis tight
        title('black: cluster [high amp]');
        subplot 312;
        plot(z1); axis tight
        title('first DSS component');
        subplot 313;
        plot(z2); axis tight
        title('last DSS component');
    else
        subplot 311;
        plot(x); hold on
        xx=x; xx(IDX{1},:)=nan;
        plot(xx,'k');
        axis tight
        title('black: cluster [high amp]');
        subplot 312;
        plot(z2); axis tight
        title('last DSS component');
        subplot 313;
        nt_sgram(z2,128,1); axis tight
        title('last DSS component');
    end
    
    if 1 
        figure(105); clf
        subplot 211;
        %sr=44100; ERBspect(z1,sr); 
        nt_sgram(z1,1024,32,[],1);
        title('first');
        subplot 212;
        %ERBspect(z2,sr);
        nt_sgram(z2,1024,32,[],1);
        title('last');
    end
    if nargout==0; clear IDX SCORE TODSS; end
    
end

function y=norm2(x,n,ind)
[I,J]=ind2sub([n,n],ind);
for k=1:size(x,1)
    a=x(k,1:n);
    b=sqrt(a(I).*a(J));
    y(k,:)=x(k,:)./b;
end

    
    
    