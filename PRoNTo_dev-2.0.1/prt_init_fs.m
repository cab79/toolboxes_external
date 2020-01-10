function [fid,PRT,tocomp] = prt_init_fs(PRT, in, mids,mask,precmask,headers)
% function to initialise the kernel data structure
% ------------------------------------------------
%
% FORMAT: Two modes are possible:
%     fid = prt_init_fs(PRT, in)
%     [fid, PRT, tocomp] = prt_init_fs(PRT, in)
%
% USAGE 1:
% -------------------------------------------------------------------------
% function will return the id of a feature set or an error if it doesn't
% exist in PRT.mat
% Input:
% ------
% in.fs_name: name for the feature set (string)
%
% Output:
% -------
% fid : is the identifier for the feature set in PRT.mat
%
% USAGE 2:
% -------------------------------------------------------------------------
% function will create the feature set in PRT.mat and overwrite it if it
% already exists.
% Input:
% ------
% in.fs_name: name for the feature set (string)
% in.fname:   name of PRT.mat
%
% in.mod(m).mod_name:  name of the modality
% in.mod(m).detrend:   type of detrending
% in.mod(m).mode:      'all_scans' or 'all_cond'
% in.mod(m).mask:	   mask used to create the feature set
% in.mod(m).param_dt:  parameters used for detrending (if any)
% in.mod(m).normalise: scale the input scans or not
% in.mod(m).matnorm:   mat file used to scale the input scans
%
% Output:
% -------
% fid : is the identifier for the model constructed in PRT.mat
%
% Populates the following fields in PRT.mat (copied from above):
%   PRT.fs(f).fs_name
%   PRT.fs(f).fas
%   PRT.fs(f).k_file
% Also computes the following fields:
%   PRT.fs(f).id_mat:       Identifier matrix (useful later)
%   PRT.fs(f).id_col_names: Columns in the id matrix
%
% Note: this function does not write PRT.mat. That should be done by the
%       calling function
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

% find index for the new feature set
fs_exists = false;
if ~(prt_checkAlphaNumUnder(in.fs_name))
    beep
    disp('Feature set name should be entered in alphanumeric format only')
    disp('Please correct')
    return
end
if isfield(PRT,'fs')
    if any(strcmpi(in.fs_name,{PRT.fs(:).fs_name}))
        fid = find(strcmpi(in.fs_name,{PRT.fs(:).fs_name}));
        fs_exists = true;
    else
        fid = length(PRT.fs)+1;
    end
else
    fid = 1;
end

% do we want to initialise the feature set?
if nargout == 1
    if ~fs_exists
        error('prt_init_fs:fsNotFoundinPRT',...
            ['Feature set ''',in.fs_name,''' not found in PRT.mat.']);
    end
else
    % initialise
    
    [pathName fileName]=fileparts(in.fname);
    if fs_exists
        warning('prt_init_fs:overwriteFsInPRT',...
            ['Feature set ''',in.fs_name,''' found in PRT.mat. Overwriting ...']);
    else
        % doesn't exist. initialise the structure
        disp(['Feature set ''',in.fs_name,''' not found in PRT.mat. Creating...'])
    end
    
    PRT.fs(fid).fs_name = in.fs_name;
    PRT.fs(fid).k_file = in.fs_name;
    PRT.fs(fid).id_col_names = {'group','subject','modality','condition','block','scan'};
    PRT.fs(fid).fas=struct('im',[],'ifa',[]);
    n_vox=0;
    n_mods=length(mids);
    for m = 1:n_mods
        PRT.fs(fid).modality(m).mod_name = in.mod(mids(m)).mod_name;
        PRT.fs(fid).modality(m).detrend  = in.mod(mids(m)).detrend;
        PRT.fs(fid).modality(m).param_dt = in.mod(mids(m)).param_dt;
        PRT.fs(fid).modality(m).mode     = in.mod(mids(m)).mode;
        %get indexes from mask specified in the data and design step
        vm = spm_vol(mask{m});
        vm = spm_read_vols(vm);
        if ~any(vm(:)>0)
            error('prt_init_fs:NoVoxelinMask',...
                ['Mask of modality ',num2str(m),' does not contain any voxel >0'])
        else
            PRT.fs(fid).modality(m).feat_idx_img = find(vm>0);
        end
        mid = mids(m);
        if m==1
            n_vox = sum(vm(:)>0);
        end
        if n_vox ~= sum(vm(:)>0)
            error('prt_init_fs:MasksNotConsistent',...
                'Masks access areas of different sizes across modalities')
        end
        %get subindexes from mask specified in the data prepare step
        if ~isempty(precmask{m})
            vm = spm_vol(precmask{m});
            vm = spm_read_vols(vm);
            if ~any(vm(:)>0)
                error('prt_init_fs:NoVoxelinMask',...
                    ['2nd level mask of modality ',num2str(m),' does not contain any voxel >0'])
            end
            [d,PRT.fs(fid).modality(m).idfeat_fas] = intersect(PRT.fs(fid).modality(m).feat_idx_img, find(vm~=0));
        else
            PRT.fs(fid).modality(m).idfeat_fas=[];
        end
        PRT.fs(fid).modality(m).normalise=struct('type',[],'scaling',[]);
    end
    
    indm = zeros(n_mods,1);
    szm = zeros(n_mods,1);
    
    % First count the total number of samples. Loops are needed since each
    % subject may have a variable number of scans
    n = 0;
    for gid = 1:length(PRT.group) % group
        for sid = 1:length(PRT.group(gid).subject);  % subject
            for m = 1:n_mods
                mid = mids(m);
                if strcmpi(in.mod(mid).mode,'all_scans');
                    n = n + size(PRT.group(gid).subject(sid).modality(mid).scans,1);
                elseif strcmpi(in.mod(mid).mode,'all_cond')
                    if ~isfield(PRT.group(gid).subject(sid).modality(mid).design,'conds')
                        error('prt_init_fs:fsIsAllCondModelisAllScans',...
                            ['''All conditions'' selected for modality ', num2str(m)...
                            ' but no design was specified. This syntax is invalid, '...
                            'Please use ''All Scans'' instead.']);
                    end
                    for cid = 1:length(PRT.group(gid).subject(sid).modality(mid).design.conds)    % condition
                        n = n + length(PRT.group(gid).subject(sid).modality(mid).design.conds(cid).scans);
                    end
                end
            end  % modality
        end  % subject
    end  % group
    PRT.fs(fid).id_mat = zeros(n,length(PRT.fs(fid).id_col_names));
    PRT.fs(fid).fas.im = zeros(n,1);
    PRT.fs(fid).fas.ifa= zeros(n,1); 
    
    % Count the total number of samples and set sample ids for the kernel
    % Set fas for the file arrays
    sample_range = 0;
    for gid = 1:length(PRT.group) % group
        for sid = 1:length(PRT.group(gid).subject);  % subject
            for m = 1:n_mods
                mid = mids(m);
                
                if strcmpi(in.mod(mid).mode,'all_scans')
                    n_vols_s  = size(PRT.group(gid).subject(sid).modality(mid).scans,1);
                    all_scans = 1:n_vols_s;
                    
                    % configure indices
                    sample_range = (1:n_vols_s)+max(sample_range);
                    PRT.fs(fid).id_mat(sample_range,1) = gid;
                    PRT.fs(fid).id_mat(sample_range,2) = sid;
                    PRT.fs(fid).id_mat(sample_range,3) = mid;
                    
                    if isfield(PRT.group(gid).subject(sid).modality(mid).design,'conds')
                        conds = PRT.group(gid).subject(sid).modality(mid).design.conds;
                        for cid = 1:length(conds)
                            scans  = PRT.group(gid).subject(sid).modality(mid).design.conds(cid).scans;
                            blocks = PRT.group(gid).subject(sid).modality(mid).design.conds(cid).blocks;
                            
                            PRT.fs(fid).id_mat(sample_range(scans),4) = cid;
                            PRT.fs(fid).id_mat(sample_range(scans),5) = blocks;
                            %PRT.fs(fid).id_mat(sample_range(scans),6) = 1:length(scans);
                        end
                        
                        PRT.fs(fid).id_mat(sample_range,6) = 1:length(all_scans);
                    else
                        scans  = 1:size(PRT.group(gid).subject(sid).modality(mid).scans,1);
                        PRT.fs(fid).id_mat(sample_range,6) = scans;
                    end
                    
                    sctoadd=(1:n_vols_s)+indm(m);
                    PRT.fs(fid).fas.ifa(sample_range)=sctoadd';
                    PRT.fs(fid).fas.im(sample_range)=mid*ones(n_vols_s,1);
                    %configure indices for the file array
                    indm(m)=n_vols_s+max(indm(m));
                elseif strcmpi(in.mod(mid).mode,'all_cond')
                    conds     = PRT.group(gid).subject(sid).modality(mid).design.conds;
                    n_vols_s  = size(PRT.group(gid).subject(sid).modality(mid).scans,1);
                    
                    % now loop over conditions
                    for cid = 1:length(conds)    % condition
                        scans     = PRT.group(gid).subject(sid).modality(mid).design.conds(cid).scans;
                        blocks    = PRT.group(gid).subject(sid).modality(mid).design.conds(cid).blocks;
                        n_vol_s_c = length(scans);
                        if n_vol_s_c==0
                            sample_range = 1+max(sample_range);
                            PRT.fs(fid).id_mat(sample_range,5) = 0;
                            PRT.fs(fid).id_mat(sample_range,6) = 0;
                        else
                            sample_range = (1:n_vol_s_c)+max(sample_range);
                            PRT.fs(fid).id_mat(sample_range,5) = blocks;
                            PRT.fs(fid).id_mat(sample_range,6) = scans;
                            %configure indices for the file array
                            sctoadd=scans+indm(m);
                            PRT.fs(fid).fas.ifa(sample_range)=sctoadd';
                            PRT.fs(fid).fas.im(sample_range)=mid*ones(n_vol_s_c,1);
                        end
                        
                        % configure indices
%                         sample_range = (1:n_vol_s_c)+max(sample_range);
                        PRT.fs(fid).id_mat(sample_range,1) = gid;
                        PRT.fs(fid).id_mat(sample_range,2) = sid;
                        PRT.fs(fid).id_mat(sample_range,3) = mid;
                        PRT.fs(fid).id_mat(sample_range,4) = cid;
%                         PRT.fs(fid).id_mat(sample_range,5) = blocks;
                        %PRT.fs(fid).id_mat(sample_range,6) = 1:length(sample_range);
%                         PRT.fs(fid).id_mat(sample_range,6) = scans;
                        
                        
                    end
                    %configure indices for the file array
                    indm(m)=n_vols_s+max(indm(m));
                end
                szm(m)=szm(m)+size(PRT.group(gid).subject(sid).modality(mid).scans,1);
            end  % modality
        end  % subject
    end  % group
    
    %initialize the file arrays if they do not exist already or if the
    %detrending parameters were modified
    if ~isfield(PRT,'fas');
        % initialise all modalities (not just those we're working on)
        for m = 1:length(PRT.masks)
            
            PRT.fas(m)=struct('mod_name',[],'dat',[],'detrend',[],'param_dt',[],'hdr',[]);
            PRT.fas(m).mod_name = PRT.masks(m).mod_name;
        end
    end
    tocomp=zeros(1,length(in.mod));
    prt_dir=fileparts(in.fname);
    for i=1:n_mods
        % check whether we need to recreate the file array
        if mids(i)>length(PRT.fas) ||...
                isempty(PRT.fas(mids(i)).dat) || exist(PRT.fas(mids(i)).dat.fname)==0 ||...
                PRT.fas(mids(i)).detrend ~= in.mod(mids(i)).detrend  || ...
                (isempty(PRT.fas(mids(i)).param_dt) && ~isempty(in.mod(mids(i)).param_dt)) || ...
                (~isempty(PRT.fas(mids(i)).param_dt) && isempty(in.mod(mids(i)).param_dt)) || ...
                ((~isempty(PRT.fas(mids(i)).param_dt) && ~isempty(in.mod(mids(i)).param_dt)) && ...
                PRT.fas(mids(i)).param_dt~=in.mod(mids(i)).param_dt)
            
            if mids(i)>length(PRT.fas) || isempty(PRT.fas(mids(i)).dat)
                disp(['File array does not exist for modality ''',...
                    char(in.mod(mids(i)).mod_name),'''. Creating...'])
            elseif PRT.fas(mids(i)).detrend ~= in.mod(mids(i)).detrend ...
                    && any(strcmpi(fieldnames(PRT.fas(mids(i)).dat),'fname')) ...
                    && exist(PRT.fas(mids(i)).dat.fname,'file')
                
                warning('prt_init_fs:overwriteFileArray',...
                    ['File array already exists for modality ''',...
                    char(in.mod(mids(i)).mod_name),''', but parameters ',...
                    'have changed. Re-creating ...']);
                
                delete(PRT.fas(mids(i)).dat.fname);
            end
            
            tocomp(mids(i))=1;
            %PRT.fas(mids(i)).mod_name = in.mod(mids(i)).mod_name;
            PRT.fas(mids(i)).detrend = in.mod(mids(i)).detrend;
            PRT.fas(mids(i)).param_dt = in.mod(mids(i)).param_dt;
            PRT.fas(mids(i)).hdr = headers{i};
            PRT.fas(mids(i)).idfeat_img = PRT.fs(fid).modality(i).feat_idx_img;                % index of voxels in the full image (nifti)
            datname=[prt_dir,filesep,'Feature_set_',char(in.mod(mids(i)).mod_name),'.dat'];
            PRT.fas(mids(i)).dat = file_array(...
                datname, ...                 % fname     - filename
                [szm(i),n_vox],...           % dim       - dimensions (default = [0 0] )
                spm_type('float32'), ...  % dtype     - datatype   (default = 'float')
                0, ...                       % offset    - offset into file (default = 0)
                1);                          % scl_slope - scalefactor (default = 1)
        else
            disp(['Using existing file array for modality ''', ...
                char(in.mod(mids(i)).mod_name),'''.'])
        end
        
        % check that the input .mat for the scaling have the right size
        if in.mod(mids(i)).normalise==2
            try
                load(in.mod(mids(i)).matnorm);
            catch
                error('prt_prepare_data:ScalingMatUnloadable',...
                    'Could not load the .mat file containing the scaling')
            end
            try
                szin=max(size(scaling));
            catch
                error('prt_prepare_data:ScalingNotinFile',...
                    'This file does not contain the "scaling" field required')
            end
            if szin~=szm(i)
                error('prt_prepare_data:Scalingdimensionwrong',...
                    'The dimension of the .mat file does not correspond to the number of scans in that modality')
            end
            PRT.fs(fid).modality(i).normalise.type=2;
            PRT.fs(fid).modality(i).normalise.scaling=reshape(scaling,1,szm(i));
        elseif in.mod(mids(i)).normalise==1
            PRT.fs(fid).modality(i).normalise.type=1;
        else
            PRT.fs(fid).modality(i).normalise.type=0;
        end
    end
    
    PRT.fs(fid).modality=rmfield(PRT.fs(fid).modality,'feat_idx_img');
end
