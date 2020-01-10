function dr=prt_comp_ranking_dist(v,w)
% Function to compute the distance between two ranking vectors, as detailed
% in Lampel and Moran, 2005 (in Information Retrieval, 8, 245-264).
% 
% INPUT : two ranking vectors of the same size
% OUTPUT: their distance
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Jessica Schrouff, 18/10/2012

if nargin<2
    error('prt_comp_ranking_dist:nargin',...
        'two ranking vectors should be entered, please correct')
end
nr=length(v);

if nr~=length(w)
    error('prt_comp_ranking_dist:sizeerror',...
        'the ranking vectors do not have the same size, please correct')
end

% Compute the distance
dr=0;
for i=1:nr
    for j=1:nr
        if (v(i)<v(j) && w(i)>w(j)) 
            tmp=1;
        else
            tmp=0;
        end
        dr=dr+tmp;
    end
end
dr=2*dr/(nr*(nr-1));

end