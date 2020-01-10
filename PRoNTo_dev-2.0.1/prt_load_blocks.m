function block = prt_load_blocks(filenames, bs, br)
% Load one or more blocks of data.
% This script is a effectively a wrapper function that for the routines
% that actually do the work (SPM nifti routines)
%
% The syntax is either:
%
% img = prt_load_blocks(filenames, block_size, block_range) just to specify
% continuous blocks of data
%
% or
%
% img = prt_load_blocks(filenames, voxel_index) to access non continuous
% blocks
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

if nargin <2 || nargin >3
    disp('Usage: img = prt_load_blocks(filenames, block_size, block_range)');
    disp('or')
    disp('Usage: img = prt_load_blocks(filenames, voxel_indexes)');
    return;
end

% read the image dimensions from the header
N  = nifti(filenames);
dm = size(N(1).dat);
if length(dm)==2, dm = [dm 1]; end % handling case of 2D image
n_vox = prod(dm(1:3));

if length(dm) == 3
    n_vol = 1;
else
    n_vol = dm(4);
end

% get the data
if nargin==3
    data_range = (br(1)-1)*bs+1:min(br(end)*bs,n_vox);
else
    data_range = bs;
end

block=zeros(length(data_range),length(N));
if n_vol==1
    for i=1:length(N)
        block(:,i) = N(i).dat(data_range);
    end
else
    for i=1:n_vol
        dat_r = N(1).dat(:,:,:,i);
        block(:,i) = dat_r(data_range);
    end
end
return