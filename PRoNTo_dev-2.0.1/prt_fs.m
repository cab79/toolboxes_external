function [outfile] = prt_fs(PRT,in)
% Function to build file arrays containing the (linearly detrended) data
% and compute a linear (dot product) kernel from them
%
% Inputs:
% -------
% in.fname:      filename for the PRT.mat (string)
% in.fs_name:    name of fs and relative path filename for the kernel matrix
%
% in.mod(m).mod_name:  name of modality to include in this kernel (string)
% in.mod(m).detrend:   detrend (scalar: 0 = none, 1 = linear)
% in.mod(m).param_dt:  parameters for the kernel detrend (e.g. DCT bases)
% in.mod(m).mode:      'all_cond' or 'all_scans' (string)
% in.mod(m).mask:      mask file used to create the kernel
% in.mod(m).normalise: 0 = none, 1 = normalise_kernel, 2 = scale modality
% in.mod(m).matnorm:   filename for scaling matrix
% in.mod(m).multroi    1 if one kernel per region required
% in.mod(m).atlasroi   name of the atlas to build one kernel per region
%
% in.flag_mm:   Perform multi-kernel learning (1) or not (0)? If yes, the
% kernel is saved as a cell vector, with one kernel per modality
%
% Outputs:
% --------
% Calls prt_init_fs to populate basic fields in PRT.fs(f)...
% Writes PRT.mat
% Writes the kernel matrix to the path indicated by in.fs_name
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand and J. Schrouff
% $Id$

% Configure some variables and get defaults
% -------------------------------------------------------------------------
prt_dir  = regexprep(in.fname,'PRT.mat', ''); % or: fileparts(fname);

% get the index of the modalities for which the user wants a kernel/data
n_mods=length(in.mod);
mids=[];
for i=1:n_mods
    if ~isempty(in.mod(i).mod_name)
        mids = [mids, i];
    end
end
n_mods=length(mids);

% Load mask(s) and resize if necessary
%force second-level masking by atlas if provided
for i=1:length(in.mod)
    if isfield(in.mod(i),'atlasroi')
        if ~isempty(in.mod(i).atlasroi) && isempty(in.mod(i).mask)
            in.mod(i).mask = in.mod(i).atlasroi;    %option only available for one modality
        end
    end
end

[mask,precmask,headers,PRT,ratl] = load_masks(PRT, prt_dir, in, mids);

% Initialize the file arrays, kernel and feature set parameters
[fid,PRT,tocomp] = prt_init_fs(PRT,in,mids,mask,precmask,headers);

in.tocomp = tocomp;
in.precmask = precmask;
in.fid = fid;

% Build the feature set and kernel
%--------------------------------------------------------------------------

Phi = [];
igd = [];
PRT.fs(fid).multkernelROI = 0; % Multiple kernels with an atlas
PRT.fs(fid).multkernel = 0;    % Multiple kernels from different modalities

if in.flag_mm   % One kernel per modality so need to treat them independently
    for i = 1:n_mods  % multiple modalitiestt
        % For each modality, get the corresponding ID mat and sample index
        idtk = PRT.fs(fid).id_mat(:,3) == mids(i);
        nimm = length(unique(PRT.fs(fid).id_mat(:,3) == mids(i)));
        
        %check that modalities have the same dimensions in terms of samples
        nim1 =length(unique(PRT.fs(fid).id_mat(:,3) == mids(1)));
        if nimm~= nim1
            error('prt_fs:MultKernMod_DifIm',...
                'Modalities should have the same number of samples to be considered for MKL')
        end
        addin.ID = PRT.fs(fid).id_mat(idtk,:);
        
        % If second-level, i.e. atlas-based kernels as well
        if isfield(in.mod(mids(i)),'multroi') ...
                && in.mod(mids(i)).multroi
            atl=spm_vol(ratl{i});
            %Initialize all fields and compute the feature sets if needed
            if any(in.tocomp)
                [PRT] = prt_fs_modality(PRT,in,1,addin);
            end
            [PRT,Phim,igd] = prt_compute_ROI_kernels(PRT,in,fid,mids(i),atl,addin,i);
            PRT.fs(fid).atlas_name = ratl;   
            kerns = Phim(igd);
        else
            [PRT,Phim] = prt_fs_modality(PRT,in,1,addin);
            [d1,idmax] = max(Phim);
            [d1,idmin] = min(Phim);
            min_max = find(idmax==idmin);
            if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
                igd = i;
                Phim = {Phim};
                kerns = Phim;
            else
                beep
                disp('No overlap between data and mask/atlas for at least one sample')
                disp(['Kernel ',num2str(i),' will be removed from further analysis'])
                kerns = [];
                nmodin = 1:length(PRT.fs(fid).modality);
                igm = setdiff(nmodin,i);
                PRT.fs(fid).modality = PRT.fs(fid).modality(igm);
            end
            PRT.fs(fid).atlas_name = {};
        end
        kerns = reshape(kerns,1,length(kerns));
        Phi = [Phi, kerns];
        PRT.fs(fid).modality(i).igood_kerns = igd;
    end
    %post-hoc: the ID mat should be the same for all modalities involved,
    %so only the first one will be saved
    indm=PRT.fs(fid).fas.im==1;
    PRT.fs(fid).id_mat=PRT.fs(fid).id_mat(indm,:);
    PRT.fs(fid).multkernel = 1;
else
    % Concatenate the modalities in samples or only one modality
    % Does the same as before without playing with the ID matrix
    % If second-level, i.e. atlas-based kernels as well
    
    % First check that all concatenated modalitlies have the same flag and
    % atlas, or no atlas
    if n_mods>1
        mult = zeros(n_mods,1);
        for i = 1:n_mods
            if isfield(in.mod(mids(i)),'multroi') ...
                    && in.mod(mids(i)).multroi
                mult(i) = 1;
                atl = ratl{i};
                if i ==1
                    atlmod = atl;
                end
                if ~strcmpi(atl,atlmod)
                    error('prt_fs:ConcatenatingWithDifferentAtlases',...
                        'Concatenated multiple modalities should have the same atlas');
                end
            end
        end
        if length(unique(mult))~=1
            error('prt_fs:ConcatenateModalitiesWithandWithoutAtlas',...
                'Modalities cannot be concatenated unless they all have no or the same atlas')
        end
    end
    
    atl=spm_vol(ratl{1});
    addin = struct();
    %Initialize all fields and compute the feature sets if needed
    if isfield(in.mod(mids(1)),'multroi') ...
            && in.mod(mids(1)).multroi
        if any(in.tocomp)
            [PRT] = prt_fs_modality(PRT,in,0,[]);
        end
        [PRT,Phim,igd] = prt_compute_ROI_kernels(PRT,in,fid,mids(1),atl,addin,1);
        PRT.fs(fid).atlas_name{1} = ratl{1};
        if isempty(igd)
            error('prt_fs:NoDataInMask',...
                'No overlap between data and mask/atlas for at least one sample, cannot create kernel')
        end
        kerns = Phim(igd); % Makig sure the dimensions are consistent
        kerns = reshape(kerns,1,length(kerns));
        Phi = kerns;
    else % Simply concatenate the modalities in samples
        [PRT,Phim] = prt_fs_modality(PRT,in,0,[]);
        PRT.fs(fid).multkernel = 0;
        PRT.fs(fid).atlas_name = {};
        [d1,idmax] = max(Phim);
        [d1,idmin] = min(Phim);
        min_max = find(idmax==idmin);
        if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
            igd = 1;
            Phi{1}=Phim;
        else
            error('prt_fs:NoDataInMask',...
                'No overlap between data and mask/atlas for at least one sample, cannot create kernel')
        end
    end
    PRT.fs(fid).modality(1).igood_kerns = igd;
end

clear Phim


% Save kernel and function output
% -------------------------------------------------------------------------
outfile = in.fname;
disp('Saving feature set to: PRT.mat.......>>')
disp(['Saving kernel to: ',in.fs_name,'.mat.......>>'])
fs_file = [prt_dir,in.fs_name];
if spm_check_version('MATLAB','7') < 0
    save(outfile,'-V6','PRT');
    save(fs_file,'-V6','Phi');
else
    save(outfile,'PRT');
    save(fs_file,'Phi');
end
disp('Done.')

%--------------------------------------------------------------------------
%------------------------- Private function -------------------------------
%--------------------------------------------------------------------------

function [mask, precmask, headers,PRT, ratl] = load_masks(PRT, prt_dir, in, mids)
% function to load the mask for each modality
% -------------------------------------------
n_mods   = length(mids);
mask     = cell(1,n_mods);
precmask = cell(1,n_mods);
headers  = cell(1,n_mods);
ratl     = cell(1,n_mods);
for m = 1:n_mods
    mid = mids(m);
    
    % get mask for the within-brain voxels (from data and design)
    ddmask = PRT.masks(mid).fname;
    try
        M = nifti(ddmask);
    catch %#ok<*CTCH>
        error('prt_fs:CouldNotLoadFile',...
            'Could not load mask file');
    end
    
    % get mask for the kernel if one was specified
    mfile = in.mod(mid).mask;
    if ~isempty(mfile) %&&  mfile ~= 0
        try
            precM = spm_vol(char(mfile));
        catch
            error('prt_fs:CouldNotLoadFile',...
                'Could not load mask file for preprocessing');
        end
    end
    
    % get atlas for the ROI based kernel if one was specified
    if isfield(in.mod(mid),'atlasroi')
        alfile = in.mod(mid).atlasroi;
        if ~isempty(alfile) %&&  mfile ~= 0
            try
                precA = spm_vol(char(alfile));
            catch
                error('prt_fs:CouldNotLoadFile',...
                    'Could not load mask file for preprocessing');
            end
        end
    else
        alfile=[];
    end
    
    % get header of the first scan of that modality
    if isfield(PRT,'fas') && mid<=length(PRT.fas) && ...
            ~isempty(PRT.fas(mid).dat)
        N = PRT.fas(mid).hdr;
    else
        N = spm_vol(PRT.group(1).subject(1).modality(mid).scans(1,:));
    end
    headers{m}=N;
    
    % compute voxel dimensions and check for equality if n_mod > 1
    if m == 1
        n_vox = prod(N.dim(1:3));
    elseif n_mods > 1 && n_vox ~= prod(N.dim(1:3))
        error('prt_fs:multipleModatlitiesVariableFeatures',...
            'Multiple modalities specified, but have variable numbers of features');
    end
    
    % resize the different masks if needed
    if N.dim(3)==1, Npdim = N.dim(1:2); else Npdim = N.dim; end % handling case of 2D images
    if any(size(M.dat(:,:,:,1)) ~= Npdim)
        warning('prt_fs:maskAndImagesDifferentDim',...
            'Mask has different dimensions to the image files. Resizing...');
        
        V2 = spm_vol(char(ddmask));
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
        mfile_new = ['updated_1stlevel_mask_m',num2str(mid)];
        movefile(fullfile(V2_pth,[rV2_fn,'.img']), ...
            fullfile(prt_dir,[mfile_new,'.img']));
        movefile(fullfile(V2_pth,[rV2_fn,'.hdr']), ...
            fullfile(prt_dir,[mfile_new,'.hdr']));
        PRT.masks(mid).fname = fullfile(prt_dir,[mfile_new,'.img']);
        mask{m} = PRT.masks(mid).fname;
    else
        mask{m} = ddmask;
    end
    if ~isempty(mfile) && any((precM.dim~= N.dim)) % && mfile ~= 0
        warning('prt_fs:maskAndImagesDifferentDim',...
            'Preprocessing mask has different dimensions to the image files. Resizing...');
        V2 = spm_vol(char(mfile));
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
        % if more than one 2nd level mask to resize
        nummask = 1;
        while exist(fullfile( ...
                prt_dir,['updated_2ndlevel_mask_m',num2str(mid),'_',...
                num2str(nummask),'.img']),'file')
            nummask = nummask+1;
        end
        mfile_new = ['updated_2ndlevel_mask_m',num2str(mid),...
            '_',num2str(nummask)];
        movefile(fullfile(V2_pth,[rV2_fn,'.img']), ...
            fullfile(prt_dir,[mfile_new,'.img']));
        movefile(fullfile(V2_pth,[rV2_fn,'.hdr']), ...
            fullfile(prt_dir,[mfile_new,'.hdr']));
        precmask{m} = fullfile(prt_dir,[mfile_new,'.img']);
    else
        precmask{m} = mfile;
    end
    if ~isempty(alfile) && any((precA.dim~= N.dim))
        warning('prt_fs:atlasAndImagesDifferentDim',...
            'Atlas has different dimensions to the image files. Resizing...');
        V2 = spm_vol(char(alfile));
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
        alfile_new = ['updated_atlas_',V2_fn];
        movefile(fullfile(V2_pth,[rV2_fn,'.img']), ...
            fullfile(prt_dir,[alfile_new,'.img']));
        movefile(fullfile(V2_pth,[rV2_fn,'.hdr']), ...
            fullfile(prt_dir,[alfile_new,'.hdr']));
        ratl{m} = fullfile(prt_dir,[alfile_new,'.img']);
    else
        ratl{m} = alfile;
    end
    clear M N precM V1 V2 mfile mfile_new
end


function [PRT,Phi,igd,nroi] = prt_compute_ROI_kernels(PRT,in,fid,mids,atl,addin,ind_mod)
% function to load the mask for each modality
% -------------------------------------------

%For each region, get the indexes of the voxels in the 2nd level mask
h=spm_read_vols(atl);
if ~isempty(PRT.fs(fid).modality(ind_mod).idfeat_fas)
    idt = PRT.fs(fid).modality(ind_mod).idfeat_fas; %indexes of voxels in second-level mask
else
    idt = 1:length(PRT.fas(mids).idfeat_img);
end
idm1 = PRT.fas(mids).idfeat_img(idt); %indexes of voxels in first level mask
interh = h(idm1);
roi = unique(interh(interh>0));
nroi = length(roi);
Phi=cell(nroi,1);
in.tocomp = zeros(1,length(in.tocomp));
%For each region, compute kernel and save the indexes in the image for
%further computation of the weights
PRT.fs(fid).modality(ind_mod).idfeat_img = cell(nroi,1);
igd = []; %indexes of non 0 kernels
for i=1:nroi
    disp ([' > Computing kernel: ', num2str(i),' of ',num2str(nroi),' ...'])
    addin.idvox_fas = idt(interh == roi(i));
    [PRT,Phim] = prt_fs_modality(PRT,in,1,addin);
    [d1,idmax] = max(Phim);
    [d1,idmin] = min(Phim);
    min_max = find(idmax==idmin);
    if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
        igd = [igd,i];
    else
        beep
        disp('No overlap between data and mask/atlas for at least one sample')
        disp(['Region ',num2str(i),' will be removed from further analysis'])
    end
    Phi{i}=Phim;
    %         idts = idt(interh == roi(i));
    PRT.fs(fid).modality(ind_mod).idfeat_img{i} = find(interh == roi(i)) ;
end
PRT.fs(fid).multkernelROI = 1;
if ~isempty(igd)
    PRT.fs(fid).modality(ind_mod).idfeat_img = PRT.fs(fid).modality(ind_mod).idfeat_img(igd);
    PRT.fs(fid).modality(ind_mod).num_ROI = roi(igd);
end


return
