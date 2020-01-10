function [IDX,TODSS,SCORE]=nt_cluster_jd2(x,dsr,flags,init)
%[IDX,todss,SCORE]=nt_cluster_jd2(x,dsr,flags,init) - cluster with joint diagonalization
%
%  IDX: indices of cluster ownership
%  TODSS: DSS matrices to emphasize both clusters
%  SCORE: scores for both clusters (smaller means better contrast)
%
%  x: data (time*channels)
%  dsr: downsample ratio for cross product series
%  flags: 'norm', 'norm2': give each dsr-sized slice the same weight
%  init: provide initial clustering
%
% See nt_bias_cluster, nt_cluster1D


if nargin<2; error('!'); end
if nargin<3 ||isempty(flags); flags=[]; end
if nargin<4; init=[]; end

if ~nargout; disp('entering nt_cluster_jd...'); end

% time series of cross products (terms of covariance matrix)
[xx,ind]=nt_xprod(x,'lower',dsr);

% option: give each slice the same weight (counters amplitude variations)
if find(strcmp(flags,'norm'))
    xx=nt_normrow(xx);
end
if find(strcmp(flags,'norm2'))
    xx=norm2(xx,size(x,2),ind);
end

% initial clustering of the time series of cross products
if isempty(init)
    [C,A,score]=nt_cluster1D(xx);
    [~,idx]=min(score);
    A=A(:,idx);
    % upsample the cluster ownership index 
    A=repmat(A',[dsr,1]);
    A=A(:);
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==1); 
    IDX{2}=find(A==2);
else
    IDX{1}=setdiff(1:size(x,1),init);
    IDX{2}=init;
end

% initial DSS to contrast clusters
c0=nt_cov(x);
c1=nt_cov(x(IDX{2},:));
[todss,pwr0,pwr1]=nt_dss0(c0,c1);
todss2=todss(:,[1 end]); % keep only first and last components
z=nt_mmat(x,todss2);

PLOT_FIG2=0;
if PLOT_FIG2
    figure(2);  clf; set(gcf, 'name','in nt_cluster_jd');
    A=zeros(size(x,1),1); A(IDX{1})=1;
    subplot 511; plot(x); title('data');
    subplot 512; plot(A,'.-'); title('initial cluster map');
    subplot 513; plot(z(:,1)); title('initial DSS1');
    subplot 514; plot(z(:,2)); title('initial DSS2');
    drawnow; %pause(1);
end

% iterate until stable
old_IDX=IDX;
for k=1:10

    [zz,ind]=nt_xprod(z,[],dsr);
    
    % option: give each slice the same weight (counters amplitude variations)
    if find(strcmp(flags,'norm'))
        zz=nt_normrow(zz);
    end
    if find(strcmp(flags,'norm2'))
        zz=norm2(zz,size(z,2),ind);
    end

    if 0 % expand power terms to improve distribution of small values
        EXPONENT=0.5;
        zz(1,1)=(zz(1,1)).^EXPONENT;
        zz(2,2)=(zz(2,2)).^EXPONENT;
    end
    
    % 1D clustering of each of the 3 time series
    [C,A,score]=nt_cluster1D(zz);
    [~,idx]=min(score); % choose best
    A=A(:,idx);
    A=repmat(A',[dsr,1]);
    A=A(:);
    A(end:size(x,1))=A(end);
    IDX{1}=find(A==1); 
    IDX{2}=find(A==2);

    if ~nargout; disp(['cluster sizes: ', num2str([numel(IDX{1}), numel(IDX{2})])]); end
    
    % DSS to contrast clusters
    c1=nt_cov(x(IDX{1},:));
    c0=c1+nt_cov(x(IDX{2},:));
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    z=nt_mmat(x,todss(:,[1 end])); % keep first and last
    
    if PLOT_FIG2
        figure(2);  
        subplot 515; plot(A,'.-'); title('final cluster map');
    end
    if all(size(old_IDX{1})==size(IDX{1})) && all(old_IDX{1}==IDX{1}); break; end
    old_IDX=IDX;
end 

% scores for each cluster
c1=nt_cov(x(IDX{1},:));
c0=c1+nt_cov(x(IDX{2},:));
[TODSS{1},pwr0,pwr1]=nt_dss0(c0,c1);
SCORE(1)=pwr1(end)/pwr0(end);

c1=nt_cov(x(IDX{2},:));
c0=c1+nt_cov(x(IDX{1},:));
[TODSS{2},pwr0,pwr1]=nt_dss0(c0,c1);
SCORE(2)=pwr1(end)/pwr0(end);

SCORE=SCORE(:);

if nargout==0;
    
    % no output, just plot
    disp(['cluster1: ',num2str(100*numel(find(A==1))/numel(A)), '%']);

    figure(100); clf
    subplot 311;
    plot(x); hold on
    xx=x;
    xx(IDX{1},:)=nan;
    plot(xx,'k');
    axis tight
    title('black: cluster B');
    
    z1=nt_mmat(x,todss(:,1));
    z2=nt_mmat(x,todss(:,end));

    figure(101); clf ;
    subplot 221;
    plot(pwr1./pwr0,'.-'); xlabel('component'); ylabel('score'); title('DSS cluster A vs all');
    subplot 222;
    wsize=min(1024,size(z1,1));
    nt_spect_plot(z1,wsize,[],[],1);
    hold on
    wsize=min(1024,size(z2,1));
    nt_spect_plot(z2,wsize,[],[],1);
    xlim([0 .5])
    nt_linecolors([],[1 3]);
    legend('first','last'); legend boxoff
    hold off

    z=nt_mmat(x,todss); 
    z=nt_normcol(z);
    subplot 223; imagescc(nt_cov(z(IDX{1},:))); title('cluster A'); 
    subplot 224; imagescc(nt_cov(z(IDX{2},:))); title('cluster B');

    
    figure(100);
    subplot 312;
    plot(z1); axis tight
    title('first DSS component')
    subplot 313;
    plot(z2); axis tight
    title('last DSS component');
    
%     figure(103); clf
%     e1=mean(z(find(A==1),:).^2);
%     e2=mean(z(find(A==2),:).^2);
%     plot([e1',e2'], '.-'); legend('cluster A', 'cluster B'); title ('power per component');
    
    
    if 0 
        figure(105); clf
        subplot 211;
        nt_sgram(z1,1024,32,[],1);
        title('first');
        subplot 212;
        nt_sgram(z2,1024,32,[],1);
        title('last');
    end
    clear IDX SCORE TODSS
    
end

function y=norm2(x,n,ind)
[I,J]=ind2sub([n,n],ind);
for k=1:size(x,1)
    a=x(k,1:n);
    b=sqrt(a(I).*a(J));
    y(k,:)=x(k,:)./b;
end

    
    
    