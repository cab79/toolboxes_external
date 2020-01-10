function imagescc(a, b, c)
%imagescc - plot image with symmetric scaling

if nargin==1
    C=a;
    m=max(abs(C(:)));
    imagesc(C,[-m-realmin,m+realmin]);
else
    X=a;
    Y=b;
    C=c;
    m=max(abs(C(:)));
    imagesc(X,Y,C,[-m-realmin,m+realmin]);
end