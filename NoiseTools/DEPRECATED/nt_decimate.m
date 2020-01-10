function y=nt_decimate(x,R)
%y=nt_decimate(x,R) - apply matlab decimate function to columns of matrix
% 

sz=size(x);
m=sz(1);
x=reshape(x,[m,prod(sz(2:end))]);

y=zeros(ceil(m/R),prod(sz(2:end)));

for k=1:size(y,2);
    y(:,k)=decimate(x(:,k),R)';
end

y=reshape(y,[size(y,1),sz(2:end)]);

