function [NW_roi,idfeatroi] = prt_build_region_weights(weight_fname,atlas_fname,build_im,comp_perm)
%
% Function to compute the weights for each region as specified by the atlas
% image (one value per region). Weights not in the atlas are comprised in an
% additional region with name 'others'.
%
% --------------------------------------------------------------------------
% INPUT : 
%   - weight_fname : name of the weights image
%   - atlas_fname  : name of the atlas image
%   - build_im     : set to 1 to build the resulting image.
%   - comp_perm    : set to 1 to compute the ranking for the permutations.
% OUTPUT: 
%   - image file with the normalized weights in each region, which can
%           be viewed in the results GUI as a weight image.
%   - .mat file containing the weights of each ROI, in %, pH, the
%           normalized weights for each region, in %, pHN and the list of
%           region names.
%__________________________________________________________________________
% Copyright (C) 2015 Machine Learning & Neuroimaging Laboratory

% Written by Jessica Schrouff, 19/09/2012


%select weight map
if nargin<1
    f=spm_select(1,[],'Select weight map',[],pwd,'.img');
    if isempty(f)
        disp('No weight map selected, aborting')
        return
    end
else
    f=weight_fname{1};
end

%select atlas
if nargin<2
    gi=spm_select(1,'image','Select atlas');
    if isempty(f)
        disp('No atlas selected, aborting')
        return
    end
else
    gi=atlas_fname;
end

%set flag to 1 if not specified
if nargin<3
    flag=1;
else
    flag=build_im;
end

% in case the ranking has to be computed for the permutations, get the
% names of the weight images for each permutation.
if nargin<4
    comp_perm=0;
end

if comp_perm
    [a,b]=fileparts(char(f));
    dirn=[a,filesep,'perm_',b];
    pth=pwd;
    if isdir(dirn)
        %get the names of the weight images
        cd(dirn)
        files=dir('*perm*.img');
        fp=char({files(:).name});
        fper={[repmat(dirn,length(files),1),repmat(filesep,length(files),1),fp]};
        fperm=char([fper;{f}]);    
        cd(pth)
    else
        disp('No folder containing the weight images for permutations')
        disp('Computing normalized weights for the provided image only')
        comp_perm=0;
        fperm=f;
    end
else
    fperm=f;
end
    

%resize atlas if needed
%--------------------------------------------------------------------------

%load images
V=spm_vol(f);
[xxx,bb]=fileparts(f);
nfo=length(V);
V1=spm_vol(gi);
dumb=V(1);

if ~any(dumb.dim == V1.dim)
    disp('Resizing atlas--------->>')
    %reslice
    fl_res = struct('mean',false,'interp',0,'which',1,'prefix','resized_');
    spm_reslice([dumb V1],fl_res);

    %build updated atlas
    [V1_pth,V1_fn,V1_ext] = spm_fileparts(V1.fname);
    rV1_fn = [fl_res.prefix,V1_fn];

    if strcmp(V1_ext,'.nii')
        V_in = spm_vol(fullfile(V1_pth,[rV1_fn,'.nii']));
        V_out = V_in; V_out.fname = fullfile(V1_pth,[rV1_fn,'.img']);
        spm_imcalc(V_in,V_out,'i1');
    end

    %put the files into the PRT directory
    mfile_new=['resized_',V1_fn];
    pp=spm_fileparts(f);
    movefile(fullfile(V1_pth,[rV1_fn,'.img']), ...
    fullfile(pp,[mfile_new,'.img']));
    movefile(fullfile(V1_pth,[rV1_fn,'.hdr']), ...
    fullfile(pp,[mfile_new,'.hdr']));
    g=spm_vol(fullfile(pp,[mfile_new,'.img']));
    h=spm_read_vols(g);
else
    h=spm_read_vols(V1);
end

%compute histogram
%--------------------------------------------------------------------------

for ii=1:size(fperm,1)
    
    %Get the volumes into matrices
    V=spm_vol(fperm(ii,:));
    w=zeros(V(1).dim(1)*V(1).dim(2)*V(1).dim(3),length(V));
    VV=spm_read_vols(V);
    nfold=size(VV,4)-1;
    if nfold == 0 %when only the average across folds was computed
        w=VV(:);
    else
        for i=1:nfold+1  %number of folds + average
            tmp=VV(:,:,:,i);
            w(:,i)=tmp(:);
        end
        clear tmp
    end    
    atlas=round(h(:)); % assuming integer numbers of regions
    atlas(isnan(w(:,1)))=NaN;
    
     %Compute the volume of the 'others' region
    if ii==1       
        N_other=length(find(atlas==0));
        P_oth=N_other/length(find(isnan(atlas)));
        disp(['Volume of the others region: ',num2str(P_oth)])
    end
    
    % Compute the weights and normalized weights
    disp(['Computing weights in each ROI for image ',num2str(ii)])
    [H HN SN idfeatroi] = prt_region_histogram(w, atlas);
    nr=size(H,1);
    
    % Put 0 if the fold has only NaNs
    indnan = [];
    for i=1:size(HN,2)
        if length(find(isnan(HN(:,i))))==size(HN,1) || ...
                length(find(HN(:,i)==0))==size(HN,1)
            HN(:,i) = 0;
            H(:,i) = 0;
            indnan = [indnan i];
        end
    end
    
    %Correct for the 'others' region (one time)
    if ii==1   
        %compute proportions as in PCA
        r_min=min(atlas);
        R=max(atlas);
        if r_min==0
            corr=1;
        else
            corr=0;
        end
    end
    
    %sum of weights in each region
    pH=H*100;
    oth_w=pH(1); %save the weight of the 'others' region
   
    %normalized sum of weights in each region
    inn= ~isnan(HN(:,1));
    shn=sum(HN(inn,:),1);
    pHN=(HN./repmat(shn,size(HN,1),1))*100;
    pHN(:,indnan)=0;
        
    
end

%save sorted H, HN and the list of corresponding ROIs for the 'true' image,
%and the ranking of the average weights for each permutation
W_roi=pH;
NW_roi=pHN/100;
SN=SN*100;
[a,b,c]=fileparts(dumb.fname);
% save(fullfile(a,['atlas_',b1,'_',b,'.mat']),'LR',...
%     'W_roi','NW_roi','dwn','drwn','erwn','SN','ER','P_oth','oth_w');

%build new image with the normalized weights and save values
%--------------------------------------------------------------------------

if flag
    disp('Building image of normalized weights')
    
    % if image exists, overwrite
    if exist(fullfile( ...
            a,['ROI_',b,c]),'file')
        disp('Image of normalized weights per region already exists, overwriting...')
    end
    
    %build image if flag
    img_name=[a,filesep,'ROI_',b,c];
    img4d = file_array(img_name,size(VV),'float32-le',0,1,0);
    for km=1:size(w,2)
        for r=r_min:R
            w(idfeatroi{r+corr},km)=pHN(r+corr,km);
        end
        img4d(:,:,:,km)=reshape(w(:,km),[size(VV,1),size(VV,2), size(VV,3),1]);
    end
    
    
    % Create weigths file
    %--------------------------------------------------------------------------
    disp('Creating image--------->>')
    No         = V(1).private;     % copy header
    No.dat     = img4d;            % change file_array
    No.descrip = 'Pronto weigths'; % description
    create(No);                    % write header
    disp('Done.')
end



% Previous functions that were saving all results in a separate .mat file

    %compute the rank of each region according to the weights
%     [dub,ih]=sort(pH,1,'descend');
%     [d1,dw]=sort(ih);
%     %ranking distance between each fold and the "average" fold
%     drw=zeros(1,nfold);
%     for ifold=1:nfold
%         drw(ifold)=prt_comp_ranking_dist(dw(:,ifold),dw(:,end));
%     end
%     %Expected value of the rank for each region
%     rw=zeros(nr,1);
%     for i=1:nr
%         for j=1:nr
%             tmp=length(find(dw(i,:)==j));
%             rw(i)=rw(i)+j*tmp;
%         end
%     end
%     rw=rw/nfold;
