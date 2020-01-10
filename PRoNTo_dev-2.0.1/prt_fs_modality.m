function [PRT,Phi] = prt_fs_modality(PRT,in, flag, addin)
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
%
% in.fid:      index of feature set to create
% in.tocomp:   vector of booleans indicating whether to build the feature set
% in.precmask: cell array containing the names of the second-level mask for
%           each modality to build
%
% flag:     set to 1 to compute one kernel per region as labelled in atlas
% addin:    additional inputs for this operation to optimize computation
%
% Outputs:
% --------
% Writes the kernel matrix to the path indicated by in.fs_name and the
% feature set in a file array if it needs to be computed
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand and J. Schrouff
% $Id: prt_fs.m 650 2013-02-20 13:36:58Z amarquan $

% Configure some variables and get defaults
% -------------------------------------------------------------------------
prt_dir  = regexprep(in.fname,'PRT.mat', ''); % or: fileparts(fname);
def=prt_get_defaults('fs');
fid = in.fid;

% get the indexes of the samples and of the features to use if flag is set
% to 1
if nargin>=3 && flag
    if ~isfield(addin,'ID')
        ID= PRT.fs(fid).id_mat;
        n_mods=length(in.mod);
        mids=[];
        for i=1:n_mods
            if ~isempty(in.mod(i).mod_name)
                mids = [mids, i];
            end
        end
    else
        ID = addin.ID;
        mids = unique(ID(:,3));
    end
    n_mods=length(mids);
    if isfield(addin,'idvox_fas')
        idvox = addin.idvox_fas;
    else
        if in.tocomp(mids(1)) % build everything
            idvox = 1:PRT.fas(mids(1)).dat.dim(2); %n_vox has to be the same for all concatenated modalities (version 1.1)
        else
            if ~isempty(PRT.fs(fid).modality(mids(1)).idfeat_fas)
                idvox = PRT.fs(fid).modality(mids(1)).idfeat_fas;
            else
                idvox = 1:PRT.fas(mids(1)).dat.dim(2);
            end
        end
    end
    n_vox = numel(idvox);
    nfa = 0;    
else
    % get the index of the modalities for which the user wants a kernel/data
    n_mods=length(in.mod);
    mids=[];
    for i=1:n_mods
        if ~isempty(in.mod(i).mod_name)
            mids = [mids, i];
        end
    end
    n_mods=length(mids);
    ID= PRT.fs(fid).id_mat;
    nfa = [];
    for m = 1:n_mods
        nfa   = [nfa, PRT.fas(mids(m)).dat.dim(1)];
        if in.tocomp(mids(m))
            n_vox = PRT.fas(mids(m)).dat.dim(2); %n_vox has to be the same for all concatenated modalities (version 1.1)
        else
            if ~isempty(PRT.fs(fid).modality(m).idfeat_fas)
                idvox = PRT.fs(fid).modality(m).idfeat_fas;
            else
                idvox = 1:PRT.fas(mids(m)).dat.dim(2);
            end
            n_vox = numel(idvox);
        end
    end
end



% -------------------------------------------------------------------------
% ---------------------Build file arrays and kernel------------------------
% -------------------------------------------------------------------------
n   = size(ID,1);
Phi = zeros(n);
% set memory limit
mem         = def.mem_limit;
block_size  = floor(mem/(8*3)/max([nfa, n])); % Block size (double = 8 bytes)
n_block     = ceil(n_vox/block_size);

bstart = 1; bend = min(block_size,n_vox);
if nargin<3 || ~flag
    h = waitbar(0,'Please wait while preparing feature set');
    step=1;
end
for b = 1:n_block
    if nargin<3 || ~flag
        disp ([' > preparing block: ', num2str(b),' of ',num2str(n_block),' ...'])
    end
    vox_range  = bstart:bend;
    block_size = length(vox_range);
    kern_vols  = zeros(block_size,n);
    for m=1:n_mods
        mid=mids(m);
        %Parameters for the masks and indexes of the voxels
        %-------------------------------------------------------------------
        % get the indices of the voxels within the file array mask (data &
        % design step)
        ind_ddmask = PRT.fas(mid).idfeat_img(vox_range);
        
        %load the mask for that modality if another one was specified
        if ~isempty(in.precmask{m})
            prec_mask = prt_load_blocks(in.precmask{m},ind_ddmask);
        else
            prec_mask = ones(block_size,1);
        end
        %indexes to access the file array
        indm = find(PRT.fs(fid).fas.im==mid);
        ifa  = PRT.fs(fid).fas.ifa(indm);
        indm = find(ID(:,3) == mid);
        %get the data from each subject of each group and save its linear
        %detrended version in a file array
        %-------------------------------------------------------------------
        if in.tocomp(mid)  %need to build the file array corresponding to that modality
            n_groups = length(PRT.group);
            sample_range=0;
            nfa=PRT.fas(mid).dat.dim(1);
            datapr=zeros(block_size,nfa);
            
            %get the data for each subject of each group
            for gid = 1:n_groups
                for sid = 1:length(PRT.group(gid).subject)
                    n_vols_s = size(PRT.group(gid).subject(sid).modality(mid).scans,1);
                    sample_range = (1:n_vols_s)+max(sample_range);
                    fname = PRT.group(gid).subject(sid).modality(mid).scans;
                    datapr(:,sample_range) = prt_load_blocks(fname,ind_ddmask);
                    %check for NaNs, in case of beta maps                  
                    [inan,jnan] = find(isnan(datapr(:,sample_range)));
                    if ~isempty(inan)
                        disp('Warning: NaNs found in loaded data')
                        disp('Consider updating 1st level mask for better performance')
                        for inn=1:length(inan)
                            datapr(inan(inn),sample_range(jnan(inn))) = 0;
                        end
                    end
                    
                    %detrend if necessary
                    if in.mod(mid).detrend ~= 0
                        if  isfield(PRT.group(gid).subject(sid).modality(mid).design,'TR')
                            TR = PRT.group(gid).subject(sid).modality(mid).design.TR;
                        else
                            try
                                TR = PRT.group(gid).subject(sid).modality(mid).TR;
                            catch
                                error('detrend:TRnotfound','No TR in data, suggesting that detrend is not necessary')
                            end
                        end
                        switch in.mod(mid).detrend
                            case 1
                                C = poly_regressor(length(sample_range), ...
                                        in.mod(mid).param_dt);
                            case 2
                                C = dct_regressor(length(sample_range), ...
                                        in.mod(mid).param_dt,TR);
                        end
                        R = eye(length(sample_range)) - C*pinv(C);
                        datapr(:,sample_range) = datapr(:,sample_range)*R';
                    end
                end
            end
            
            if b==1
                % Write the detrended data into the file array .dat
                namedat=['Feature_set_',char(in.mod(mid).mod_name),'.dat'];
                fpd_clean(m) = fopen(fullfile(prt_dir,namedat), 'w','ieee-le'); %#ok<AGROW> % 'a' append
                fwrite(fpd_clean(m), datapr', 'float32',0,'ieee-le');
            else
                % Append the data in file .dat
                fwrite(fpd_clean(m), datapr', 'float32',0,'ieee-le');
            end
            
            % get the data to build the kernel
            kern_vols(:,indm) = datapr(:,ifa).* ...
                repmat(prec_mask~=0,1,length(ifa));
            % if a scaling was entered, apply it now
            if ~isempty(PRT.fs(fid).modality(m).normalise.scaling)
                kern_vols(:,indm) = kern_vols(:,indm)./ ...
                    repmat(PRT.fs(fid).modality(m).normalise.scaling,block_size,1);
            end
            clear datapr
        else
            kern_vols(:,indm) = (PRT.fas(mid).dat(ifa,idvox(vox_range)))'; % IF CRASHES HERE, RE-START MATLAB
            % if a scaling was entered, apply it now
            if ~isempty(PRT.fs(fid).modality(m).normalise.scaling)
                kern_vols(:,indm) = kern_vols(:,indm)./ ...
                    repmat(PRT.fs(fid).modality(m).normalise.scaling,block_size,1);
            end
            
        end
        if nargin<3 || ~flag
            waitbar(step/ (n_block*n_mods),h);
            step=step+1;
        end
    end
    if size(kern_vols,2)>6e3
        % Slower way of estimating kernel but using less memory.
        % size limit setup such that no more than ~1Gb of mem is required:
        % 1Gb/3(nr of matrices)/8(double)= ~40e6 -> sqrt -> 6e3 element
        for ic=1:size(kern_vols,2)
            Phi(:,ic) = Phi(:,ic) + kern_vols' * kern_vols(:,ic);
        end
    else
        Phi = Phi + (kern_vols' * kern_vols);
    end
    bstart = bend+1; bend = min(bstart+block_size-1,n_vox);
    clear block_mask kern_vols
end
if nargin<3 || ~flag
    close(h)
end

% closing feature file(s)
if exist('fpd_clean','var')
    for ii=1:numel(fpd_clean)
        fclose(fpd_clean(ii));
    end
end


% -------------------------------------------------------------------------
% Private functions
% -------------------------------------------------------------------------
function c = poly_regressor(n,order)
% n:     length of the series
% order: the order of polynomial function to fit the tend

basis = repmat([1:n]',[1 order]);
o = repmat([1:order],[n 1]);
c = [ones(n,1) basis.^o];
return

function c=dct_regressor(n,cut_off,TR)
% n:       length of the series
% cut_off: the cut off perioed in second (1/ cut off frequency)
% TR:      TR

if cut_off<0
    error('cut off cannot be negative')
end

T = n*TR;
order = floor((T/cut_off)*2)+1;
c = spm_dctmtx(n,order);
return
