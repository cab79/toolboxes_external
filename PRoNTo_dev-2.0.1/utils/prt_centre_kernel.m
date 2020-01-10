function [C,Cs,Css] = prt_centre_kernel(K, Ks, Kss)

% This function centres the kernel matrix, respecting the independence of
% training and test partitions. See Shawe-Taylor and Cristianini for
% background on this approach.
%
% Shawe-Taylor, J. and Cristianini, N. (2004). Kernel methods for Pattern
% analysis. Cambridge University Press.
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by D. Hardoon, A. Marquand and J. Mourao-Miranda
% Id: $

l = size(K,1);
j = ones(l,1);
C = K - (j*j'*K)/l - (K*j*j')/l + ((j'*K*j)*j*j')/(l^2);

if( nargin > 1 )
    tk =  (1/l)*sum(K,1); % (1 x l) 
    tl = ones(size(Ks,1),1); % (n x 1)
    Cs = Ks - ( tl * tk); % ( n x l )
    tk2 = (1/(size(Ks,2)))*sum(Cs,2); % ( n x 1 )   
    Cs = Cs - (tk2 * j'); % ( n x l )
    
    % Two equivalent ways to achieve the same thing
    %Cs = Ks - repmat(sum(K),size(Ks,1),1)/l - repmat(sum(Ks,2),1,size(Ks,2))/size(Ks,2) + repmat(j'*K*j,size(Ks,1),size(Ks,2))/(l^2);
    %Cs = Ks - (tl*sum(K))/l - (sum(Ks,2)*j')/size(Ks,2) + ((j'*K*j)*tl*j')/(l^2);
    
    if nargin > 2
        ttj = ones(size(Kss,1),1);
        Css = Kss - (sum(Ks,2)*ttj')/l - (ttj*sum(Ks,2)')/l + ((j'*K*j)*ttj*ttj')/(l^2);
    end
end 

end