function [nmask] = prt_utils_update_mask(files,mask,dir)

% Script to update the first level mask to be entered in the Data and
% Design step of a PRoNTo analysis. The mask is first resized to match the
% dimensions of provided images. Any NaN value found in any of the images
% will then be removed from the mask.
%
% Modalities/runs that will be concatenated in samples should be selected
% at the same time to ensure common features>
%
% Inputs (optional):
% - files: character array of the full file names to load
% - mask : full name of the mask to update (e.g. SPMnoeyes in /masks)
% - dir  : directory where to save the updated mask
%
% Outputs: the updated mask, saved in the specified directory.
%--------------------------------------------------------------------------
% Written by J. Schrouff, 2015, as a preprocessing step for PRoNTo.

if nargin<1
    files=spm_select([1 Inf],'image','Select files to be analyzed in PRoNTo');
end

if nargin<2
    mask = spm_select(1,'image','Select maks image to be updated');
end

if nargin<3
    dir = uigetdir();
end


try
    M = nifti(mask);
catch %#ok<*CTCH>
    error('prt_utils_update_mask:CouldNotLoadFile',...
        'Could not load mask file');
end

% First, resize mask to 1st image dimensions
% -------------------------------------------------------------------------
try 
    N = spm_vol(files(1,:));
catch
    error('prt_utils_update_mask:CouldNotLoadFile',...
        'Could not load first image file');
end

if N.dim(3)==1, Npdim = N.dim(1:2); else Npdim = N.dim; end % handling case of 2D images
if any(size(M.dat(:,:,:,1)) ~= Npdim)
    warning('prt_utils_update_mask:maskAndImagesDifferentDim',...
        'Mask has different dimensions to the image files. Resizing...');
    
    V2 = spm_vol(char(mask));
    % reslicing V2
    fl_res = struct('mean',false,'interp',0,'which',1,'prefix','tmp_');
    spm_reslice([N V2],fl_res)
    % now renaming the file
    [V2_pth,V2_fn,V2_ext] = spm_fileparts(V2.fname);
    rV2_fn = [fl_res.prefix,V2_fn];
    if strcmp(V2_ext,'.nii')
        % turn .nii into .img/.hdr image!
        V_in = spm_vol(fullfile(V2_pth,[rV2_fn,'.nii']));
        V_out = V_in; V_out.fname = fullfile(V2_pth,[rV2_fn,'.img']);
        spm_imcalc(V_in,V_out,'i1');
    end
    mfile_new = ['resized_',V2_fn];
    movefile(fullfile(V2_pth,[rV2_fn,'.img']), ...
        fullfile(dir,[mfile_new,'.img']));
    movefile(fullfile(V2_pth,[rV2_fn,'.hdr']), ...
        fullfile(dir,[mfile_new,'.hdr']));
    mask = fullfile(dir,[mfile_new,'.img']);
end
V2 = spm_vol(mask);
M = spm_read_vols(V2);

% Second, loop over all files to get the voxels with NaN values
% -------------------------------------------------------------------------
uNan = [];
for i=1:size(files,1)
    try
        N = spm_vol(files(i,:));
    catch
        error('prt_utils_update_mask:CouldNotLoadFile',...
            ['Could not load image file ',num2str(i)]);
    end
    voxval = spm_read_vols(N);
    voxval = voxval(:);
    vNan = find(isnan(voxval));
    if ~isempty(vNan)
        uNan = union(uNan,vNan);
    end
end

if ~isempty(uNan)
    M(uNan) = 0;
end

[V2_pth,V2_fn,V2_ext] = spm_fileparts(V2.fname);
mfile_new = ['updated_mask_',V2_fn];
nmask = V2;
nmask.fname = fullfile(dir,[mfile_new,V2_ext]);
nmask = spm_write_vol(nmask,M);

