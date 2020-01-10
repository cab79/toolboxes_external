function [H HN SN idfeatroi] = prt_region_histogram(beta, atlas)

% 
% Atlas-based region histograms.
% 
% Based on the following function:
% 
% L1-HISTOGRAM
% (c) Luca Baldassarre
% CS Dept, UCL, London, UK
% 8th May 2012
% l.baldassarre@cs.ucl.ac.uk
% baldazen@gmail.com
%
% Atlas-based region histograms.
%
% For each column of beta, l1_histogram(beta, atlas) computes the relative amount of
% the l1_norm that is contained in each region defined by the atlas.
% Atlas is a nx1 vector, where n = size(beta,1), such that atlas(i) is the
% region to which voxel i belongs.
%
% H = l1_histogram(beta, atlas) only computes the standard histogram
%
% [H HN] = l1_histogram(beta, atlas) also computes the normalized
% histogram, where each bin is normalized by the region's volume (i.e. the
% number of voxels it contains).
%
% [H HN sorted_regions] = l1_histogram(beta, atlas) return the list of
% regions, sorted in descending order according to the normalized
% histogram.
%__________________________________________________________________________
% Copyright (C) 2015 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$


% Number of vectors
m = size(beta, 2);
% Number of regions
R = max(atlas);
% Initial region index
r_min = min(atlas);
% Add an offset to account for matlab indexing (it starts from 1)
if r_min == 0
   correction = 1;
else
   correction = 0;
end

H = zeros(R,m);
S = zeros(R,m);
for km = 1:m
    disp(['Fold ',num2str(km)])
   % Compute relative frequencies for each region
   for r = r_min:R
       tmp = sum(abs(beta(atlas == r,km)));
       H(r+correction,km) = tmp;

      %compute the proportions of positive versus negative weights
      S(r+correction,km) = length(find(beta(atlas == r,km)>0));
   end
end

%% COMPUTE NORMALIZED HISTOGRAMS AND FULL INTERSECTION
if nargout > 1
   % Compute volumes according to atlas
   volume = zeros(R,1);
   idfeatroi = cell(R,1);
   for r = r_min:R
      volume(r+correction) = sum(atlas == r);
      idfeatroi{r+correction} = find(atlas==r);
   end
   HN = H./repmat(volume,1,m);
   SN = S./repmat(volume,1,m);
end
