function data = prt_cfg_design
% Data & design configuration file
% This builds the PRT.mat data and design structure.
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M.J.Rosa
% $Id$

% ---------------------------------------------------------------------
% covar Covariates
% ---------------------------------------------------------------------
covar         = cfg_files;
covar.tag     = 'covar';
covar.name    = 'Covariates';
covar.help    = {['Select a .mat file containing '...
    'your covariates (i.e. any other data/information '...
    'you would like to include in your design). This file '...
    'should contain a variable ''R'' with a matrix of '...
    'covariates. On covariate per image is expected.']};
covar.val{1}  = {''};
covar.filter  = 'mat';
covar.ufilter = '.*';
covar.num     = [0 1];

% ---------------------------------------------------------------------
% rt_subj One per subject/scans
% ---------------------------------------------------------------------
rt_subj         = cfg_entry;
rt_subj.tag     = 'rt_subj';
rt_subj.name    = 'Regression targets (per scans)';
rt_subj.help    = {['Enter one regression target per scan. '...
    'or enter the name of a variable. '...
    ' This variable should be a vector '...
    '[Nscans x 1], where Nscans is the number of '...
    'scans/images.']};
rt_subj.strtype = 'e';
rt_subj.val     = {[]};
rt_subj.num     = [Inf 0];

% ---------------------------------------------------------------------
% regtrial One per trial
% ---------------------------------------------------------------------
rt_trial         = cfg_entry;
rt_trial.tag     = 'rt_trial';
rt_trial.name    = 'Regression targets (trials)';
rt_trial.help    = {['Enter one regression target per trial. '...
    'This vector should have the following dimensions: '...
    '[Ntrials x 1], where Ntrials is the number of trials.']
    };
rt_trial.strtype = 'e';
rt_trial.val     = {[]};
rt_trial.num     = [Inf 0];

% ---------------------------------------------------------------------
% TR Interscan interval
% ---------------------------------------------------------------------
TR         = cfg_entry;
TR.tag     = 'TR';
TR.name    = 'Interscan interval';
TR.help    = {'Specify interscan interval (TR). The units should be seconds.'};
TR.strtype = 'e';
TR.num     = [Inf 1];

% ---------------------------------------------------------------------
% units Units for design
% ---------------------------------------------------------------------
unit         = cfg_menu;
unit.tag     = 'unit';
unit.name    = 'Units for design';
unit.help    = {['The onsets of events or blocks can be specified in '...
    'either scans or seconds.']};
unit.labels  = {
    'Scans'
    'Seconds'
    }';
unit.values  = {0 1};
unit.val     = {1};

% ---------------------------------------------------------------------
% review Review
% ---------------------------------------------------------------------
review         = cfg_menu;
review.tag     = 'review';
review.name    = 'Review';
review.help    = {['Choose ''Yes'' if you would like to review your '...
    'data and design in a separate window.']};
review.labels  = {
    'No'
    'Yes'
    }';
review.values  = {0 1};
review.val     = {0};

% ---------------------------------------------------------------------
% hrfover HRF Overlap
% ---------------------------------------------------------------------
hrfover         = cfg_entry;
hrfover.tag     = 'hrfover';
hrfover.name    = 'HRF overlap';
hrfover.help    = {['If using fMRI data please specify the width of the '...
    'hemodynamic response function (HRF). This will be '...
    'used to calculate the overlap between events. '...
    'Leave as 0 for other modalities (other than fMRI).']};
hrfover.strtype = 'e';
hrfover.num     = [1 1];
hrfover.def     = @(val)prt_get_defaults('datad.hrfw', val{:});

% ---------------------------------------------------------------------
% hrfdel HRF Delay
% ---------------------------------------------------------------------
hrfdel         = cfg_entry;
hrfdel.tag     = 'hrfdel';
hrfdel.name    = 'HRF delay';
hrfdel.help    = {['If using fMRI data please specify the delay of the '...
    'hemodynamic response function (HRF). This will be '...
    'used to calculate the overlap between events. Leave '...
    'as 0 for other modalities (other than fMRI).']};
hrfdel.strtype = 'e';
hrfdel.num     = [1 1];
hrfdel.def     = @(val)prt_get_defaults('datad.hrfd', val{:});

% ---------------------------------------------------------------------
% fmri_des fMRI design specific parameters
% ---------------------------------------------------------------------
fmri_des      = cfg_branch;
fmri_des.tag  = 'fmri_des';
fmri_des.name = 'fMRI_Des';
fmri_des.val  = {hrfover, hrfdel};
fmri_des.help = {'fMRI design specific parameters, HRF overlap and delay.'};

% ---------------------------------------------------------------------
% mod_name Name
% ---------------------------------------------------------------------
mod_name         = cfg_entry;
mod_name.tag     = 'mod_name';
mod_name.name    = 'Name';
mod_name.help    = {['Name of modality. Example: ''BOLD''. The names '...
    'should be consistent accross subjects/groups '...
    'and the same names specified in the masks.']};
mod_name.strtype = 's';
mod_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {['Select scans (images) for this modality. They must '...
    'all have the same image dimensions, orientation, '...
    'voxel size etc.']};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];

% ---------------------------------------------------------------------
% subjects Subjects
% ---------------------------------------------------------------------
subjects         = cfg_files;
subjects.tag     = 'subjects';
subjects.name    = 'Files';
subjects.help    = {['Select scans (images) for this modality. They must '...
    'all have the same image dimensions, orientation, '...
    'voxel size etc.']};
subjects.filter  = 'image';
subjects.num     = [0 Inf];

% ---------------------------------------------------------------------
% modality Modality
% ---------------------------------------------------------------------
modality      = cfg_branch;
modality.tag  = 'modality';
modality.name = 'Modality';
modality.val  = {mod_name, subjects, rt_subj, covar };
modality.help = {'Specify modality, such as name and data.'};

% ---------------------------------------------------------------------
% images Scans
% ---------------------------------------------------------------------
images        = cfg_repeat;
images.tag    = 'images';
images.name   = 'Scans';
images.values = {modality };
images.help   = {['Depending on the type of data at hand, you may have many images (scans) '...
    'per subject, such as a fMRI time series, or you may have many '...
    'subjects with only one or a small number of images (scans) per subject, '...
    'such as PET images. ',...
    'Select this option if you have many subjects per modality to spatially ',...
    'normalise, but there is one or a small number of scans for '...
    'each subject. This is a faster option with less information to specify '...
    'than the ''select by subjects'' option. Both options create the same '...
    '''PRT.mat'' but ''select by scans'' is optimised for modalities '...
    'with no design.']};

% ---------------------------------------------------------------------
% fmask File name
% ---------------------------------------------------------------------
fmask        = cfg_files;
fmask.tag    = 'fmask';
fmask.name   = 'File';
fmask.filter = 'image';
fmask.ufilter = '.*';
fmask.num    = [1 1];
fmask.help   = {['Select one first-level mask (image) for each modality. ',...
    'This mask is used to optimise the prepare data step. ',...
    'In ''specify model'' there is an option to enter a ',...
    'second-level mask, which might be used to select only ',...
    'a few areas of the brain for subsequent analyses.']};
% ---------------------------------------------------------------------
% mask Modality
% ---------------------------------------------------------------------
mask         = cfg_branch;
mask.tag     = 'mask';
mask.name    = 'Modality';
mask.help    = {['Specify name of modality and file for each mask. ',...
    'The name should be consistent with the names chosen ',...
    'for the modalities (subjects/scans).']};
mask.val     = {mod_name, fmask };

% ---------------------------------------------------------------------
% masks Masks
% ---------------------------------------------------------------------
masks         = cfg_repeat;
masks.tag     = 'masks';
masks.name    = 'Masks';
masks.help    = {['Select first-level (pre-processing) mask for each ',...
    'modality. The name of the modalities should be the same ',...
    'as the ones entered for subjects/scans.']};
masks.num     = [1 Inf];
masks.values  = {mask };

% ---------------------------------------------------------------------
% load_SPM Load SPM.mat
% ---------------------------------------------------------------------
load_SPM         = cfg_files;
load_SPM.tag     = 'load_SPM';
load_SPM.name    = 'Load SPM.mat';
load_SPM.help    = {['Load design from SPM.mat (if you have previously '...
    'specified the experimental design with SPM).']};
load_SPM.filter  = '^SPM\.mat$';
load_SPM.num     = [1 1];

% ---------------------------------------------------------------------
% cond_name Name
% ---------------------------------------------------------------------
cond_name         = cfg_entry;
cond_name.tag     = 'cond_name';
cond_name.name    = 'Name';
cond_name.help    = {'Name of condition (alphanumeric strings only).'};
cond_name.strtype = 's';
cond_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% onsets Onsets
% ---------------------------------------------------------------------
onsets         = cfg_entry;
onsets.tag     = 'onsets';
onsets.name    = 'Onsets';
onsets.help    = {'Specify a vector of onset times for this condition type. '};
onsets.strtype = 'e';
onsets.num     = [Inf 1];

% ---------------------------------------------------------------------
% durations Durations
% ---------------------------------------------------------------------
durations         = cfg_entry;
durations.tag     = 'durations';
durations.name    = 'Durations';
durations.help    = {['Specify the event durations. Epoch and '...
    'event-related responses are modeled in exactly '...
    'the same way but by specifying their different '...
    'durations.  Events are specified with a duration '...
    'of 0.  If you enter a single number for the '...
    'durations it will be assumed that all trials '...
    'conform to this duration. If you have multiple '...
    'different durations, then the number must match '...
    'the number of onset times.']};
durations.strtype = 'e';
durations.num     = [Inf 1];

% ---------------------------------------------------------------------
% conds Condition
% ---------------------------------------------------------------------
conds         = cfg_branch;
conds.tag     = 'conds';
conds.name    = 'Condition';
conds.help    = {'Specify condition: name, onsets and duration.'};
conds.val     = {cond_name, onsets, durations};

% ---------------------------------------------------------------------
% conditions Conditions
% ---------------------------------------------------------------------
conditions         = cfg_repeat;
conditions.tag     = 'conditions';
conditions.name    = 'Conditions';
conditions.help    = {['Specify conditions. You are allowed to combine '...
    'both event- and epoch-related responses in '...
    'the same model and/or regressor. Any number of '...
    'condition (event or epoch) types can be '...
    'specified.  Epoch and event-related responses '...
    'are modeled in exactly the same way by '...
    'specifying their onsets [in terms of onset '...
    'times] and their durations.  Events are specified '...
    'with a duration of 0.  If you enter a single '...
    'number for the durations it will be assumed that '...
    'all trials conform to this duration.For factorial '...
    'designs, one can later associate these experimental '...
    'conditions with the appropriate levels of experimental '...
    'factors.']};
conditions.values  = {conds};

% ---------------------------------------------------------------------
% multi_conds Multiple conditions
% ---------------------------------------------------------------------
multi_conds         = cfg_files;
multi_conds.tag     = 'multi_conds';
multi_conds.name    = 'Multiple conditions';
multi_conds.val{1}  = {''};
multi_conds.help    = {
    'Select the *.mat file containing details of your multiple experimental conditions. '
    ''
    'If you have multiple conditions then entering the details a condition at a time is very inefficient. This option can be used to load all the required information in one go. You will first need to create a *.mat file containing the relevant information. '
    ''
    'This *.mat file must include the following cell arrays (each 1 x n): names, onsets and durations. eg. names=cell(1,5), onsets=cell(1,5), durations=cell(1,5), then names{2}=''SSent-DSpeak'', onsets{2}=[3 5 19 222], durations{2}=[0 0 0 0], contain the required details of the second condition. These cell arrays may be made available by your stimulus delivery program, eg. COGENT. The duration vectors can contain a single entry if the durations are identical for all events.'
    ''
    'Time and Parametric effects can also be included. For time modulation include a cell array (1 x n) called tmod. It should have a have a single number in each cell. Unused cells may contain either a 0 or be left empty. The number specifies the order of time modulation from 0 = No Time Modulation to 6 = 6th Order Time Modulation. eg. tmod{3} = 1, modulates the 3rd condition by a linear time effect.'
    ''
    'For parametric modulation include a structure array, which is up to 1 x n in size, called pmod. n must be less than or equal to the number of cells in the names/onsets/durations cell arrays. The structure array pmod must have the fields: name, param and poly.  Each of these fields is in turn a cell array to allow the inclusion of one or more parametric effects per column of the design. The field name must be a cell array containing strings. The field param is a cell array containing a vector of parameters. Remember each parameter must be the same length as its corresponding onsets vector. The field poly is a cell array (for consistency) with each cell containing a single number specifying the order of the polynomial expansion from 1 to 6.'
    ''
    'Note that each condition is assigned its corresponding entry in the structure array (condition 1 parametric modulators are in pmod(1), condition 2 parametric modulators are in pmod(2), etc. Within a condition multiple parametric modulators are accessed via each fields cell arrays. So for condition 1, parametric modulator 1 would be defined in  pmod(1).name{1}, pmod(1).param{1}, and pmod(1).poly{1}. A second parametric modulator for condition 1 would be defined as pmod(1).name{2}, pmod(1).param{2} and pmod(1).poly{2}. If there was also a parametric modulator for condition 2, then remember the first modulator for that condition is in cell array 1: pmod(2).name{1}, pmod(2).param{1}, and pmod(2).poly{1}. If some, but not all conditions are parametrically modulated, then the non-modulated indices in the pmod structure can be left blank. For example, if conditions 1 and 3 but not condition 2 are modulated, then specify pmod(1) and pmod(3). Similarly, if conditions 1 and 2 are modulated but there are 3 conditions overall, it is only necessary for pmod to be a 1 x 2 structure array.'
    ''
    'EXAMPLE:'
    'Make an empty pmod structure: '
    '  pmod = struct(''name'',{''''},''param'',{},''poly'',{});'
    'Specify one parametric regressor for the first condition: '
    '  pmod(1).name{1}  = ''regressor1'';'
    '  pmod(1).param{1} = [1 2 4 5 6];'
    '  pmod(1).poly{1}  = 1;'
    'Specify 2 parametric regressors for the second condition: '
    '  pmod(2).name{1}  = ''regressor2-1'';'
    '  pmod(2).param{1} = [1 3 5 7]; '
    '  pmod(2).poly{1}  = 1;'
    '  pmod(2).name{2}  = ''regressor2-2'';'
    '  pmod(2).param{2} = [2 4 6 8 10];'
    '  pmod(2).poly{2}  = 1;'
    ''
    'The parametric modulator should be mean corrected if appropriate. Unused structure entries should have all fields left empty.'
    }';
multi_conds.filter  = 'mat';
multi_conds.ufilter = '.*';
multi_conds.num     = [0 1];

% ---------------------------------------------------------------------
% new_design Specify design
% ---------------------------------------------------------------------
new_design         = cfg_branch;
new_design.tag     = 'new_design';
new_design.name    = 'Specify design';
new_design.help    = {'Specify design: scans (data), onsets and durations.'};
new_design.val     = {unit conditions multi_conds}; %covar for covar per trial (v3)

% ---------------------------------------------------------------------
% no_design No design
% ---------------------------------------------------------------------
no_design         = cfg_const;
no_design.tag     = 'no_design';
no_design.name    = 'No design';
no_design.val     = {0};
no_design.help    = {['Do not specify design. This option can be used '...
    'for modalities (e.g. structural scans) that do not '...
    'have an experimental design.']};

% ---------------------------------------------------------------------
% design Data & Design
% ---------------------------------------------------------------------
design        = cfg_choice;
design.tag    = 'design';
design.name   = 'Data & Design';
design.help   = {'Specify data and design.'};
design.values = {load_SPM, new_design, no_design };
design.val    = {load_SPM };

% ---------------------------------------------------------------------
% subject Modality
% ---------------------------------------------------------------------
subject      = cfg_branch;
subject.tag  = 'subject';
subject.name = 'Modality';
subject.val  = {mod_name, TR, scans, design };
subject.help = {'Add new modality.'};

% ---------------------------------------------------------------------
% gr_name Name
% ---------------------------------------------------------------------
gr_name         = cfg_entry;
gr_name.tag     = 'gr_name';
gr_name.name    = 'Name';
gr_name.help    = {'Name of the group. Example: ''Controls''.'};
gr_name.strtype = 's';
gr_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% ind_subj Subject
% ---------------------------------------------------------------------
ind_subj         = cfg_repeat;
ind_subj.tag     = 'ind_subj';
ind_subj.name    = 'Subject';
ind_subj.help    = {'Add new modality for this subject.'};
ind_subj.values  = {subject };

% ---------------------------------------------------------------------
% subjs Subjects
% ---------------------------------------------------------------------
subjs         = cfg_repeat;
subjs.tag     = 'subjs';
subjs.name    = 'Subjects';
subjs.help    = {'Add subjects/scans.'};
subjs.values  = {ind_subj };

% ---------------------------------------------------------------------
% select Select by
% ---------------------------------------------------------------------
select        = cfg_choice;
select.tag    = 'select';
select.name   = 'Select by';
select.values = {subjs, images};
select.help   = {...
    ['Depending on the type of data at hand, you may have many images (scans) '...
    'per subject, such as a fMRI time series, or you may have many '...
    'subjects with only one or a small number of images (scans) per subject '...
    ', such as PET images. If you have many scans per subject select the '...
    'option ''subjects''. If you have one scan for many subjects select '...
    'the option ''scans''.']};

% ---------------------------------------------------------------------
% group Group
% ---------------------------------------------------------------------
group         = cfg_branch;
group.tag     = 'group';
group.name    = 'Group';
group.help    = {'Specify data and design for the group.'};
group.val     = {gr_name, select };

% ---------------------------------------------------------------------
% groups Groups
% ---------------------------------------------------------------------
groups         = cfg_repeat;
groups.tag     = 'groups';
groups.name    = 'Groups';
groups.help    = {['Add data and design for one group. Click ''new'' '...
    'or ''repeat'' to add another group.']};
groups.num     = [1 Inf];
groups.values  = {group };

% ---------------------------------------------------------------------
% dir_name Directory
% ---------------------------------------------------------------------
dir_name         = cfg_files;
dir_name.tag     = 'dir_name';
dir_name.name    = 'Directory';
dir_name.help    = {['Select a directory where the PRT.mat file '...
    'containing the specified design and data matrix '...
    'will be written.']};
dir_name.filter  = 'dir';
dir_name.ufilter = '.*';
dir_name.num     = [1 1];

% ---------------------------------------------------------------------
% data Data & Design
% ---------------------------------------------------------------------
data        = cfg_exbranch;
data.tag    = 'data';
data.name   = 'Data & Design';
data.val    = {dir_name groups masks fmri_des review};
data.help   = {'Specify the data and design for each group (minimum one group).'};
data.prog   = @prt_run_design;
data.vout   = @vout_data;
data.check  = @check_data;

%------------------------------------------------------------------------
%% Output function
%------------------------------------------------------------------------
function cdep = vout_data(job)
% Specifies the output from this modules, i.e. the filename of the mat file

cdep(1)            = cfg_dep;
cdep(1).sname      = 'PRT.mat file';
cdep(1).src_output = substruct('.','files');
cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
ngroup = length(job.group);
for g = 1:ngroup
    cdep(1+g)            = cfg_dep; %#ok<*AGROW>
    cdep(1+g).sname      = sprintf('Group#%d name',g);
    cdep(1+g).src_output = substruct('.',sprintf('gr_name%d',g));
    cdep(1+g).tgt_spec   = cfg_findspec({{'strtype','s'}});
end
nmod = length(job.mask);
for m = 1:nmod
    cdep(1+ngroup+m)            = cfg_dep;
    cdep(1+ngroup+m).sname      = sprintf('Mod#%d name',m);
    cdep(1+ngroup+m).src_output = substruct('.',sprintf('mod_name%d',m));
    cdep(1+ngroup+m).tgt_spec   = cfg_findspec({{'strtype','s'}});
end

%------------------------------------------------------------------------
%% Checking functions
%------------------------------------------------------------------------
function t = check_data(data)
% Checking that the data are consistent.
t   = {};
if isfield(data,'group')
    if isfield(data,'mask')
        % Checking the names of modalities in masks and groups
        % get modality names from all subjects, if any
        Ngroup = numel(data.group);
        subj_mod = {}; gr_mod = cell(1,Ngroup); Ngr_mod = zeros(1,Ngroup);
        for gg = 1:numel(data.group)
            c_smod = 0;
            if isfield(data.group(gg).select,'subject')
                Nsubj = numel(data.group(gg).select.subject);
                for ss = 1:Nsubj
                    if ~isempty(data.group(gg).select.subject{ss})
                        for mm = 1:numel(data.group(gg).select.subject{ss})
                            c_smod = c_smod+1;
                            subj_mod{c_smod,gg} = ...
                                data.group(gg).select.subject{ss}(mm).mod_name;
                        end
                    end
                end
                gr_mod{:,gg} = unique(subj_mod(1:Nsubj,gg));%[afm]
            elseif isfield(data.group(gg).select,'modality')
                Nmod = numel(data.group(gg).select.modality);
                for mm = 1:Nmod
                    c_smod = c_smod+1;
                    subj_mod{c_smod,gg} = data.group(gg).select.modality(mm).mod_name;
                end
                gr_mod{:,gg} = unique(subj_mod(:,gg));
            end
            Ngr_mod(gg) = numel(gr_mod{:,gg});
        end
        % get modality names from all masks, if any
        Nmask = numel(data.mask);
        mask_mod = cell(1,Nmask);
        for mk = 1:Nmask
            mask_mod{mk} = data.mask(mk).mod_name;
        end
        % make sure all groups have same #modalities
        if Ngroup>1
            if any(diff(Ngr_mod))
                t = {'Different number of modalities across groups!'};
                warndlg(t,'Group modality name');
                return
            end
        end
        % NOTE: should we check if all groups have the same modalities??
        
        % make sure masks' name fit those of the subjects'
        % and no duplicate mask names
        if numel(unique(mask_mod))~=Nmask
            t = {'Masks'' name not unique!'};
            warndlg(t,'Masks modality name');
            return
        end
        ok = zeros(Nmask,1);
        for mk=1:Nmask
            if any(strcmp(mask_mod{mk},subj_mod(:)))
                ok(mk) = 1;
            end
        end
        if ~all(ok)
            t{1} = 'Some masks name don''t match any subj/img modality!';
            warndlg(t,'Masks and data modality name');
            return
        end
        smods = subj_mod(:);
        smods = smods(~cellfun(@isempty,smods));
        usubj_mod = unique(smods);
        Nsubj_mod = numel(usubj_mod);
        ok = zeros(Nsubj_mod,1);
        for sm = 1:Nsubj_mod
            if any(strcmp(usubj_mod{sm},mask_mod))
                ok(sm) = 1;
            end
        end
        if ~all(ok)
            t{1} = 'Some subj/img modality! don''t match any masks name';
            warndlg(t,'Masks and data modality name');
            return
        end
        
    end
end
return;
%-------------------------------------------------------------------------