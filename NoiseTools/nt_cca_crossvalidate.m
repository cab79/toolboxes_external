function [AA,BB,RR]=nt_cca_crossvalidate(xx,yy,lags)
%[AA,BB,RR]=nt_cca_crossvalidate(xx,yy,lags) - CCA with cross-validation
%
%  AA, BB: cell arrays of transform matrices
%  RR: r scores (2D)
%
%  xx,yy: cell arrays of column matrices
%  lags: array of lags to apply to y relative to x (can be negative)

if nargin<3; lags=[0]; end
if nargin<2; error('!'); end
if ~iscell(xx) || ~iscell(yy); error('!'); end
if length(xx) ~= length (yy); error('!'); end
if size(xx{1},1) ~= size(yy{1},1); error('!'); end

%%
% calculate covariance matrices
nTrials=length(xx);
n=size(xx{1},2)+size(yy{1},2);
C=zeros(n,n,length(lags),nTrials);
disp('Calculate all covariances...');
nt_whoss;
for iTrial=1:nTrials
    C(:,:,:,iTrial)=nt_cov_lags(xx{iTrial}, yy{iTrial},lags);
end

%%
% calculate leave-one-out CCAs
disp('Calculate CCAs...');
for iOut=1:nTrials
    CC=sum(C(:,:,:,setdiff(1:nTrials,iOut)),4); % covariance of all trials except iOut
    [A,B,R]=nt_cca([],[],[],CC,size(xx{1},2));  % corresponding CCA
    AA{iOut}=A;
    BB{iOut}=B;
end
clear C CC

%%
% calculate leave-one-out correlation coefficients
clear r rr;
disp('Calculate cross-correlations...');
for iOut=1:nTrials
    A=AA{iOut};
    B=BB{iOut};
    for iShift=1:length(lags)
        [x,y]=nt_relshift(xx{iOut},yy{iOut},lags(iShift));
        a=A(:,:,iShift);
        b=B(:,:,iShift);
        r(:,iShift)=diag( nt_normcol(x*a)' * nt_normcol(y*b )) / size(x,1); 
    end
    RR(:,:,iOut)=r;
end
disp('done');

%%
% If no output arguments, plot something informative

if nargout==0
    figure(1); clf;
    if length(lags)>1; 
        plot(mean(RR,3)'); title('correlation for each CC'); xlabel('lag'); ylabel('correlation');
    else
        plot(squeeze(mean(RR,3))); title ('correlation for each CC'); xlabel('CC'); ylabel('correlation');
    end
    figure(2); clf;
    for k=1:min(4,size(RR,2))
        subplot(2,2,k);
        [~,idx]=max(mean(RR(k,:,:),3));
        [x,y]=nt_relshift(xx{1},yy{1},lags(idx));
        plot([x*A(:,k,idx), y*B(:,k,idx)]);
        disp(corr(nt_normcol([x*A(:,k,idx), y*B(:,k,idx)])));
        title(['CC ',num2str(k)]); xlabel('sample'); 
    end
end


