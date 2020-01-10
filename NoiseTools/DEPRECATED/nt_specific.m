function y=nt_specific(x,w,thresh)
%y=nt_specific(x) - isolate channel-specific activity
%
%  y: channel-specific component
%  
%  x: data (2D or 3D)
%  w: weights
%  thresh: threshold to discard PCs (default: 10^-6)

if nargin<3||isempty(thresh); thresh=0; end
if nargin<2; w=[]; end

if ~isempty(w); error('weights not yet implemented'); end

if ndims(x)>2;
    [m,n,o]=size(x);
    x=nt_specific(nt_unfold(x),w,thresh);
    x=nt_fold(x,m);
    return
end

[m,n]=size(x);
cc=x'*x;

y=zeros(size(x));
for k=1:n
    idx=[1:k-1,k+1:n];
    [topcs,eigenvalues]=nt_pcarot(cc(idx,idx)); % PCA to orthogonalize the other channels
    topcs=topcs(:,find(eigenvalues/max(eigenvalues)>thresh));
    b=x(:,idx)*topcs(:,:);
    b=nt_normcol(b); % this could be optimized
    c=x(:,k)'*b/m; % projection matrix
    %y(:,k)=x(:,k)-b*c'; % remove projection
    y(:,k)=b*c'; % remove projection
end



% function x=nt_specific(x,w,thresh)
% %y=nt_specific(x) - isolate channel-specific activity
% %
% %  y: channel-specific component
% %  
% %  x: data (2D or 3D)
% %  w: weights
% %  thresh: threshold to discard PCs (default: 10^-6)
% 
% if nargin<3||isempty(thresh); thresh=0; end
% if nargin<2; w=[]; end
% 
% if ~isempty(w); error('weights not yet implemented'); end
% 
% if ndims(x)>2;
%     [m,n,o]=size(x);
%     x=nt_specific(nt_unfold(x),w,thresh);
%     x=nt_fold(x,m);
%     return
% end
% [m,n]=size(x);
% cc=x'*x;
% 
% y=zeros(size(x));
% for k=1:n
%     idx=[1:k-1,k+1:n];
%     [topcs,eigenvalues]=nt_pcarot(cc(idx,idx)); % PCA to orthogonalize the other channels
%     topcs=topcs(:,1:60);
%     %topcs=topcs(:,find(eigenvalues/max(eigenvalues)>thresh));
%     b=x(:,idx)*topcs(:,:);
%     b=nt_normcol(b); % this could be optimized
%     c=x(:,k)'*b/m; % projection matrix
%     x(:,k)=x(:,k)-b*c'; % remove projection
%     cc(idx,k)=x(:,k)'*x(:,idx); % update covariance matrix
%     cc(k,:)=cc(:,k)';
% end
