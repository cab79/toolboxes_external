function x=nt_enforce_stationarity(x,DSR,thresh);
% y=nt_enforce_stationarity(x,DSR,thresh) - locally project out non-stationary components
%
% y: cleaned data
%
% x: data to clean (time * channels)
% DSR: size of chunks over which to cluster
% thresh: variance ratio threshold

if nargin<3 || isempty(thresh); thresh=10; end
if nargin<2 || isempty(DSR); error('!'); end

nx=size(x,1);
x=nt_unfold(x);


score=inf;

figure(10); clf; set(gcf, 'name', 'nt_enforce_stationarity');
subplot 413; plot(x); title('original')

while 1
    
    % find clusters with maximally different covariance
    [A]=nt_cluster_jd((nt_pca(x)),DSR);
    figure(10);
    subplot 411; plot(A,'.-'); title('cluster mask'); drawnow
    
    % DSS to enhance the smaller cluster
    if numel(find(A==1))<numel(find(A==2));
        idx=find(A==1);
        [todss]=nt_dss0(nt_cov(x),nt_cov(x(find(A==1),:)));
    else
        idx=find(A==2);
        %todss=fliplr(todss);
        [todss]=nt_dss0(nt_cov(x),nt_cov(x(find(A==2),:)));
    end
    
    z=nt_mmat(x,todss);
    subplot 412; plot(z(:,1)); title('component to remove');
    score=mean(z(idx,1).^2)/mean(z(:,1).^2);
    disp(['power ratio score: ', num2str(score)]);
    
    if score<thresh
        break
    end
    
    x(idx,:)=nt_tsr_nodemean(x(idx,:),z(idx,1));
    subplot 414; plot(x); title('clean');
    drawnow
    
    %pause;
    %x=nt_demean(x);
end
x=nt_fold(x,nx);
if nargout==0; 
    clear x
end

    
    