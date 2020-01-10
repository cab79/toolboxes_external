function [SLres,Pout,XYZmm,D] = crc_parSL(Pin,opt)

% Function to hack through some PRT structure (from PRoNTo) in order to
% proceed with a "search light" approach.
%
% The idea is simple:
% 1/ proceed with your whole brain PRoNTo classification
% 2/ use the full prt structure as a template for "search light" and
%    loop over the voxels in the *1st level mask*
%       - create a small spherical volume
%       - re-estimate the kernel from the subset of voxels
%       - launch the model estimation
%       - collect tha accuracies
% 3/ save accuracies (full, balanced & per class) into images
%
% FORMAT
% [SLres,Pout] = crc_SL(Pprt,opt)
%
% INPUT:
% - Pin     Filename (with path) of PRT.mat file to use
% - opt     some options (more could be added...)
%       .R       radius in mm of the spherical searchlight [10 by def.]
%       .i_model index of PRoNTo model to use [1 by def.]
%       .loadF   load all features or not [true by def.]
%       .savImg  save results in image format or not [true by def.]
%       .permStat assess accuracy through permuation [false by def.]
%
% OUTPUT:
% - SLres   Searchlight results structure array. There is 1 structure per
%           voxel in the mask + 1 (last one) with the original whole mask
%           results
% - Pout    Filenames of generated stuff
%
%--------------------------------------------------------------------------
% NOTE:
% - PRoNTo must be initialized before using the script as it relies on
%   PRoNTo's machinery.
% - other clique format could be used, for example: cubes in images or
%   planes/lines for other types of data (connectivity matrix, ERP's, etc.)
% - no inference at the moment but this could be added...
%__________________________________________________________________________
% Copyright (C) 2015 Cyclotron Research Centre

% Written by Christophe Phillips
% University of Liège, Belgium

% Initialize PRoNTo if needed
if ~exist('prt_checkAlphaNumUnder','file')
    prt('startup','nogui')
end

% Deal with inputs
opt_def = struct( ...
    'R',10, ...
    'i_model',1, ...
    'loadF',true, ...
    'savImg',false, ...
    'plot',false, ...
    'permStat',false);
if nargin<2
    opt = opt_def;
else
    opt = prt_check_flag(opt_def,opt);
end
R = opt.R; % search light radius in mm
i_model = opt.i_model; % Model index to use
loadF = opt.loadF;
savImg = opt.savImg;
permStat = opt.permStat;
ploton = opt.plot;

if nargin<1 || isempty(Pin)
    Pprt = spm_select(1,'mat','Select the PRT.mat file');
else
    Pprt = Pin;
end
[pth,nam,ext,num] = spm_fileparts(Pprt); %#ok<*NASGU,*ASGLU>
load(Pprt)

% fix paths in PRT in case analysis directory was moved
[~,fname,ext] = fileparts(PRT.fas.dat.fname);
PRT.fas.dat.fname = fullfile(pth,[fname ext]);
save(Pprt,'PRT');

PRTorig = PRT; %#ok<*NODEF>
Pmsk = PRT.masks.fname ; % should be the updated 1st level mask
Vmsk = spm_vol(Pmsk);

%-Get space details
M         = Vmsk.mat;                          %-voxels to mm matrix
iM        = inv(M);                             %-mm to voxels matrix
DIM       = Vmsk.dim;                          %-image dimensions
[x,y,z]   = ndgrid(1:DIM(1),1:DIM(2),1:DIM(3));
XYZ       = [x(:)';y(:)';z(:)']; clear x y z    %-voxel coordinates {vx}
XYZmm     = M(1:3,:)*[XYZ; ones(1,size(XYZ,2))];%-voxel coordinates {mm}
XYZmm_cpy = XYZmm;                              %-copy without masking

% List of mask voxels>1
lVx = PRT.fas.idfeat_img;

%-Search volume (from mask)
MM     = spm_read_vols(Vmsk);
MM(isnan(MM)) = 0; % turn NaN's into 0's
MM     = logical(MM);
XYZmm  = XYZmm(:,MM(:));
XYZ    = XYZ(:,MM(:));
c            = round(DIM(:)/2);

% sanity check
ll = find(MM(:));
% any(lVx-ll) % It should be 0

switch opt.searchtype
    case {'spacetime','3Dspace'}
    %-Searchlight options (clique definition)
    searchopt.def = 'sphere';
    searchopt.spec = R;
    nSx = numel(lVx);
    
    case 'time'
    %-Searchlight options (clique definition)
    searchopt.def = 'box';
    searchopt.spec = R;
    %[~,XYZtime_firstind,XYZtime_allind] = unique(XYZ(3,:));
    XYZtimeind = find(XYZ(1,:)==c(1) & XYZ(2,:)==c(2)); % get mid-points
    XYZ = XYZ(:,XYZtimeind);
    XYZmm = XYZmm(:,XYZtimeind);
    %XYZ(:,1) = c(1); 
    %XYZ(:,2) = c(2); 
    %XYZmm = XYZmm(:,XYZtime_firstind);
    nSx = length(XYZtimeind);
end

%-Build local clique
xY      = searchopt;
xY.xyz  = [NaN NaN NaN];
xY.rej  = {'cluster','mask'};
xY      = spm_ROI(xY);

%-Build local clique, in voxels
xY.xyz       = M(1:3,:) * [c;1];
[xY, clique] = spm_ROI(xY,XYZmm_cpy);
clique       = round(iM(1:3,:) * [clique;ones(1,size(clique,2))]);
clique       = bsxfun(@minus, clique, c);
dc           = (max(clique,[],2) - min(clique,[],2) + 1)'; % size of clique radii, in voxels

%% Get the input ready before looping
PRTw = rmfield(PRTorig,'model'); % PRT before model estimation, still need to adjust feature selection!
PRTw.model = rmfield(PRTorig.model,'output');

if loadF
    % Load all features in memory
    fs_whole = PRTw.fas.dat(:,:);
else
    % or use filearray -> slower but less memory hungry!
    fs_whole = PRTw.fas.dat;
end

% find feature set index: i_fs
fs_names = cellstr(char(PRTorig.fs(:).fs_name));
i_fs = find(strcmp(PRTorig.model(i_model).input.fs.fs_name,fs_names));
if isempty(i_fs)
    error('prtSL:fs_model','No matching feature set for model #%d.\n',i_model)
end

%% Loops over all voxels from the 1st level mask & collect accuracies
% initialier SL results structure array and include at N+1 the full mask results
tmp_stats = PRTorig.model(i_model).output.stats;
if ~permStat && isfield(PRTorig.model(i_model).output.stats,'permutation')
    tmp_stats = rmfield(tmp_stats,'permutation');
end
SLres(nSx+1) = tmp_stats;
kern_file = fullfile(pth,PRTorig.fs(i_fs).k_file);

checkp = gcp('nocreate')
if opt.parallel
    if isempty(checkp)
        myPool = parpool
    end
    parforArg = Inf;
else
    parforArg = 0;
end

tic
parfor (isx = 1:nSx,parforArg)
%parfor ivx=1:nVx
    [i_SLres,igd] = process_ROI(opt.searchtype,PRTw,lVx,isx,clique,XYZ,i_model,i_fs,fs_whole,kern_file,Pprt,DIM,permStat);
    if igd
        SLres(isx) = i_SLres;
    end
end
toc

D.Vmsk = Vmsk;
D.pth = pth;
D.PRTw = PRTw;
D.DIM = DIM;
D.nVx = length(XYZ);
D.i_model = i_model;
D.R = R;
D.lVx = lVx;
fn_SLres = sprintf('SL_R%d_results.mat',R);
P_SLresults = spm_file(Pprt,'filename',fn_SLres);
save(P_SLresults,'SLres','XYZmm','D')

Pout{1} = P_SLresults;

end

%% SUBFUNCTIONS
%%=========================================================================

%% Process a single searchlight-PRT model
function [i_SLres,igd] = process_ROI(searchtype,PRT,lVx,ivx,clique,XYZ,i_model,i_fs,fs_whole,kern_file,Pprt,DIM,permStat)

%i_vx = lVx(ivx);
% 1/ update the list of voxels for 2nd level mask in PRT
xyz_cl = bsxfun(@plus,XYZ(:,ivx),clique);
% remove stuff outside image volume
l_2rm = logical(sum(bsxfun(@lt, xyz_cl, ones(3,1))));
xyz_cl(:,l_2rm) = [];
l_2rm = logical(sum(bsxfun(@gt, xyz_cl, DIM')));
xyz_cl(:,l_2rm) = [];
% create list
Lcl = sub2ind(DIM,xyz_cl(1,:),xyz_cl(2,:),xyz_cl(3,:));
Lres = find(any(bsxfun(@eq, lVx', Lcl')));
% List of voxels to use
PRT.fs(i_fs).modality.idfeat_fas = lVx(Lres);

% 2/ Rebuild kernel & save it in new file
datSL = fs_whole(:,Lres);
Phim = datSL*datSL';
[~,idmax] = max(Phim);
[~,idmin] = min(Phim);
min_max = find(idmax==idmin);
if isempty(min_max) || unique(Phim(:,min_max))~=0 %Kernel does not contain a whole line of zeros
    igd = true; % good voxel -> proceed with model estimation
    Phi{1} = Phim;
    kern_f_ivx = [kern_file,'_',turn_num2char(ivx)];
    save(kern_f_ivx,'Phi');
    PRT.fs(i_fs).k_file = spm_file(kern_f_ivx,'basename');
    %     save(kern_file,'Phi');
    
    % 3/ launch estimate
    %     PRT = prt_model(PRT,mod_w);
    in.fname      = Pprt;
    in.model_name = PRT.model(i_model).model_name;
    in.savePRT = false;
    
    [~,PRT] = prt_cv_model(PRT, in);
    %load(PRTfname)
    
    % 4/ collect results
    i_SLres = PRT.model(i_model).output.stats;
    
    % TO DO
    % 5/ permuation testing to get a p-value
    % -> only do this if it's worth (check the accuracy) !
    if permStat
        %if isfield(job,'perm_test') % to ensure back compatibility with older batch
        %    if isfield(job.perm_test,'perm_t')
        %        if isfield(job.perm_test.perm_t,'flag_sw') %keep compatibility
        %            flag = job.perm_test.perm_t.flag_sw;
        %        else
        %            flag = 0;
        %        end
                flag = 0; % save weights
                fname = Pprt;
                [~,PRT]=prt_permutation(PRT, permStat, i_model, ...
                    spm_str_manip(fname,'h'),flag);
        %    end
        %end
        i_SLres.permutation = PRT.model(i_model).output.stats.permutation;
    end
        
    % 6/ remove temporary kernel file
    delete([kern_f_ivx,'.mat']);
else
    igd = false;
    fprintf('\n Skipping voxel %d.',ivx);
    i_SLres = [];
end

end

%% turn index into chars for file numbering
function i_char = turn_num2char(i_num,nMxZ)
% Function turning an index into a char, with zero padding
if nargin<2
    % Upper limit on number to convert, by default:
    % 1e10 - 1 = 999999999, i.e. 9 digits
    nMxZ = 9;
    Zpad = '000000000';
else
    nMxZ = round(nMxZ);
    Zpad = repmat('0',1,nMxZ);
end

i_num = round(i_num);
i_char = num2str(i_num);
i_char = [Zpad(1:(nMxZ-length(i_char))) , i_char];

end
