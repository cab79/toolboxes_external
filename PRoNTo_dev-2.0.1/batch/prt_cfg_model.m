function model = prt_cfg_model
% Data & design configuration file
% This configures the kernel construction for each modality.
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Andre Marquand
% $Id$

def = prt_get_defaults;

% ---------------------------------------------------------------------
% filename Filename(s) of data
% ---------------------------------------------------------------------
infile        = cfg_files;
infile.tag    = 'infile';
infile.name   = 'Load PRT.mat';
infile.ufilter = 'PRT.mat';
infile.num    = [1 1];
infile.help   = {'Select data/design structure file (PRT.mat).'};

% ---------------------------------------------------------------------
% model_name Name
% ---------------------------------------------------------------------
model_name         = cfg_entry;
model_name.tag     = 'model_name';
model_name.name    = 'Model name';
model_name.help    = {'Name for model'};
model_name.strtype = 's';
model_name.num     = [1 Inf];
 
% ---------------------------------------------------------------------
% use_kernel Use Kernels
% ---------------------------------------------------------------------
use_kernel         = cfg_menu;
use_kernel.tag     = 'use_kernel';
use_kernel.name    = 'Use kernels';
use_kernel.help    = {...
    ['Are the data for this model in the form of kernels/basis functions? ', ...
     'If ''No'' is selected, it is assumed the data are in the form of ',...
     'feature matrices']};
use_kernel.labels  = {
               'Yes'
               'No'
}';
use_kernel.values  = {1 0};
use_kernel.val     = {1};

% ---------------------------------------------------------------------
% all_features All features
% ---------------------------------------------------------------------
all_features         = cfg_const;
all_features.tag     = 'all_features';
all_features.name    = 'All Features';
all_features.val     = {1};
all_features.help    = {...
    'Include all features from all modalities in this feature set'};

% ---------------------------------------------------------------------
% fs_name Feature set name
% ---------------------------------------------------------------------
fs_name         = cfg_entry;
fs_name.tag     = 'fs_name';
fs_name.name    = 'Name';
fs_name.help    = {'Name of a feature set. Must match design specification'};
fs_name.strtype = 's';
fs_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% mod_name Modality name
% ---------------------------------------------------------------------
mod_name         = cfg_entry;
mod_name.tag     = 'mod_name';
mod_name.name    = 'Modality name';
mod_name.help    = {'Name of modality. Example: ''BOLD''. Must match design specification'};
mod_name.strtype = 's';
mod_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% fset Feature set
% ---------------------------------------------------------------------
% fset         = cfg_branch;
% fset.tag     = 'fset';
% fset.name    = 'Feature set';
% fset.help    = {'Feature set to include in this model'};
% fset.val     = {fs_name};
            
% ---------------------------------------------------------------------
% fsets Feature sets
% ---------------------------------------------------------------------
fsets         = cfg_entry;
fsets.tag     = 'fsets';
fsets.name    = 'Feature sets';
fsets.help    = {['Enter the name of a feature set to include in this model. ',...
                  'This can be kernel or a feature matrix. ', ...
                  ]};
fsets.num     = [1 Inf];
%fsets.values  = {fsets};
fsets.strtype  = 's';

% ---------------------------------------------------------------------
% gr_name Group name
% ---------------------------------------------------------------------
gr_name         = cfg_entry;
gr_name.tag     = 'gr_name';
gr_name.name    = 'Group name';
gr_name.help    = {'Name of the group to include. Must exist in PRT.mat'};
gr_name.strtype = 's';
gr_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% all_cond All conditions
% ---------------------------------------------------------------------
all_cond         = cfg_const;
all_cond.tag     = 'all_cond';
all_cond.name    = 'All Conditions';
all_cond.val     = {1};
all_cond.help    = {'Include all conditions in this model'};

% ---------------------------------------------------------------------
% cond_name Condition name
% ---------------------------------------------------------------------
cond_name         = cfg_entry;
cond_name.tag     = 'cond_name';
cond_name.name    = 'Name';
cond_name.help    = {'Name of condition to include.'};
cond_name.strtype = 's';
cond_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% conds Conditions
% ---------------------------------------------------------------------
conds         = cfg_branch;
conds.tag     = 'conds';
conds.name    = 'Condition';
conds.help    = {'Specify condition:.'};
conds.val     = {cond_name};

% ---------------------------------------------------------------------
% sel_cond Select conditions
% ---------------------------------------------------------------------
sel_cond         = cfg_repeat;
sel_cond.tag     = 'sel_cond';
sel_cond.name    = 'Specify Conditions';
sel_cond.help    = {'Specify the name of conditions to be included '};
sel_cond.values  = {conds};

% ---------------------------------------------------------------------
% all_scans All scans
% ---------------------------------------------------------------------
all_scans         = cfg_const;
all_scans.tag     = 'all_scans';
all_scans.name    = 'All scans';
all_scans.val     = {1};
all_scans.help    = {['No design specified. This option can be used '...
    'for modalities (e.g. structural scans) that do not '...
    'have an experimental design or for an fMRI design',...
    'where you want to include all scans in the timeseries']};

% ---------------------------------------------------------------------
% conditions Conditions
% ---------------------------------------------------------------------
conditions        = cfg_choice;
conditions.tag    = 'conditions';
conditions.name   = 'Conditions / Scans';
conditions.values = {sel_cond, all_cond, all_scans};
conditions.help   = {...
    ['Which task conditions do you want to include? '...
    'Select conditions: select specific conditions from the timeseries. ', ...
    'All conditions: include all conditions extracted from the timeseries. ', ...
    'All scans: include all scans for each subject. This may be used for ', ...
    'modalities with only one scan per subject (e.g. PET), ', ...
    'if you want to include all scans from an fMRI timeseries (assumes you ',...
    'have not already detrended the timeseries and extracted task components)']};

% ---------------------------------------------------------------------
% modality Modality
% ---------------------------------------------------------------------
% modality      = cfg_branch;
% modality.tag  = 'modality';
% modality.name = 'Modality';
% modality.val  = {mod_name conditions};
% modality.help = {'Specify modality, such as name and data.'};

% ---------------------------------------------------------------------
% modalities Modalities
% ---------------------------------------------------------------------
% modalities         = cfg_repeat;
% modalities.tag     = 'modalities';
% modalities.name    = 'Modalities';
% modalities.help    = {...
%     ['Add modalities. Note that if multiple modalities are entered here, ',...
%      'they will be added to the feature set as additional samples. ',...
%      'In other words, modalities are concatenated along the sample ',...
%      'dimension, not the feature dimension. For example, this is ',...
%      'appropriate for accommodating multiple fMRI runs from identical ',...
%      'subjects.']};
% modalities.num     = [1 Inf];
% modalities.values  = {modality};

% ---------------------------------------------------------------------
% subj_num Subjects selected (per group)
% ---------------------------------------------------------------------
subj_nums         = cfg_entry;
subj_nums.tag     = 'subj_nums';
subj_nums.name    = 'Subjects';
subj_nums.help    = {
    ['Subject numbers to be included in this class. Note that individual ',...
     'numbers (e.g. 1), or a range of numbers ',...
     '(e.g. 3:5) can be entered'] };
subj_nums.strtype = 'e';
subj_nums.num     = [Inf 1];

% ---------------------------------------------------------------------
% group Group
% ---------------------------------------------------------------------
group         = cfg_branch;
group.tag     = 'group';
group.name    = 'Group';
group.help    = {'Specify data and design for the group.'};
group.val     = {gr_name, subj_nums, conditions};

% ---------------------------------------------------------------------
% class_name Class name
% ---------------------------------------------------------------------
class_name         = cfg_entry;
class_name.tag     = 'class_name';
class_name.name    = 'Name';
class_name.help    = {'Name for this class, e.g. ''controls'' '};
class_name.strtype = 's';
class_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% groups Groups
% ---------------------------------------------------------------------
groups         = cfg_repeat;
groups.tag     = 'groups';
groups.name    = 'Groups';
groups.help    = {['Add one group to this class. Click ''new'' '...
                    'or ''repeat'' to add another group.']};
groups.num     = [1 Inf];
groups.values  = {group};

% ---------------------------------------------------------------------
% class Class
% ---------------------------------------------------------------------
class         = cfg_branch;
class.tag     = 'class';
class.name    = 'Class';
class.help    = {...
    ['Specify which groups, modalities, subjects and conditions should ',...
     'be included in this class']};
class.val     = {class_name, groups};

% ---------------------------------------------------------------------
% reg_targets Regression Targets
% ---------------------------------------------------------------------
reg_targets         = cfg_entry;
reg_targets.tag     = 'reg_targets';
reg_targets.name    = 'Regression targets';
reg_targets.help    = {['Specify continuous valued target variables']};
reg_targets.strtype = 'e';
reg_targets.num     = [Inf 1];

% % ---------------------------------------------------------------------
% % mod_name Modality name
% % ---------------------------------------------------------------------
% mod_name2         = cfg_entry;
% mod_name2.tag     = 'mod_name2';
% mod_name2.name    = 'Modality name';
% mod_name2.help    = {'Name of modality. We only allow one modality for regression model per group at this moment' ...
%     'Example: ''BOLD''. Must match design specification'};
% mod_name2.strtype = 's';
% mod_name2.num     = [1 Inf];

% ---------------------------------------------------------------------
% reg_group Regression group
% ---------------------------------------------------------------------
reg_group         = cfg_branch;
reg_group.tag     = 'reg_group';
reg_group.name    = 'Group';
reg_group.help    = {'Specify data and design for the group.'};
%reg_group.val     = {gr_name, subj_nums, conditions, reg_targets};
%reg_group.val     = {gr_name, subj_nums,mod_name2 };
reg_group.val     = {gr_name, subj_nums};

% ---------------------------------------------------------------------
% k_args Define k for partioning
% ---------------------------------------------------------------------
k_args         = cfg_entry;
k_args.tag     = 'k_args';
k_args.name    = 'k';
k_args.help    = {['Number of folds/partitions for CV. To create a 50%-50%,' ...
    'choose k as 2. Please note that there can be more partitions than'...
    ' specified when leaving subjects per group out. Also note that '...
    'leaving more than 50% of the data out is not permitted.']};
k_args.strtype = 'e';
k_args.num     = [1 1];

% ---------------------------------------------------------------------
% cv_loso Leave-one-subject-out
% ---------------------------------------------------------------------
cv_loso         = cfg_const;
cv_loso.tag     = 'cv_loso';
cv_loso.name    = 'Leave one subject out';
cv_loso.val     = {1};
cv_loso.help    = {'Leave a single subject out each cross-validation iteration'};

% ---------------------------------------------------------------------
% cv_lkso K-folds CV on Subjects
% ---------------------------------------------------------------------
cv_lkso         = cfg_branch;
cv_lkso.tag     = 'cv_lkso';
cv_lkso.name    = 'k-folds CV on subjects';
cv_lkso.val     = {k_args};
cv_lkso.help    = {'k-partitioning of subjects at each cross-validation iteration'};

% ---------------------------------------------------------------------
% cv_losgo Leave-one-subject-per-group-out
% ---------------------------------------------------------------------
cv_losgo         = cfg_const;
cv_losgo.tag     = 'cv_losgo';
cv_losgo.name    = 'Leave one subject per group out';
cv_losgo.val     = {1};
cv_losgo.help    = {...
    ['Leave out a single subject from each group at a time. ', ...
     'Appropriate for repeated measures or paired samples designs.']};
 
% ---------------------------------------------------------------------
% cv_lksgo K_folds CV on Subjects per Group
% ---------------------------------------------------------------------
cv_lksgo         = cfg_branch;
cv_lksgo.tag     = 'cv_lksgo';
cv_lksgo.name    = 'k-folds CV on subjects per group';
cv_lksgo.val     = {k_args};
cv_lksgo.help    = {...
    ['K-partitioning of subjects from each group at a time. ', ...
     'Appropriate for repeated measures or paired samples designs.']};
 
% ---------------------------------------------------------------------
% cv_lobo Leave-one-block-out
% ---------------------------------------------------------------------
cv_lobo         = cfg_const;
cv_lobo.tag     = 'cv_lobo';
cv_lobo.name    = 'Leave one block out';
cv_lobo.val     = {1};
cv_lobo.help    = {...
    ['Leave out a single block or event from each subject each iteration. ', ...
     'Appropriate for single subject designs.']};

% ---------------------------------------------------------------------
% cv_lkbo K-fold CV on blocks
% ---------------------------------------------------------------------
cv_lkbo         = cfg_branch;
cv_lkbo.tag     = 'cv_lkbo';
cv_lkbo.name    = 'k-folds CV on blocks';
cv_lkbo.val     = {k_args};
cv_lkbo.help    = {...
    ['k-partitioning on blocks or events from each subject each iteration. ', ...
     'Appropriate for single subject designs.']};
 
% ---------------------------------------------------------------------
% cv_loro Leave--one-run-per-subject-out (leave one modality out per
% subject)
% ---------------------------------------------------------------------
cv_loro         = cfg_const;
cv_loro.tag     = 'cv_loro';
cv_loro.name    = 'Leave one run/session out';
cv_loro.val     = {1};
cv_loro.help    = {...
    ['Leave out a single run (modality) from each subject each iteration. ', ...
     'Appropriate for single subject designs with multiple runs/sessions.']};
   
% ---------------------------------------------------------------------
% cv_custom Feature set mask
% ---------------------------------------------------------------------
cv_custom        = cfg_files;
cv_custom.tag    = 'cv_custom';
cv_custom.name   = 'Custom';
cv_custom.filter = 'mat';
cv_custom.ufilter = '.*';
cv_custom.num    = [1 1];
cv_custom.help   = {...
    'Load a cross-validation matrix comprising a CV variable'};

% ---------------------------------------------------------------------
% cv_type Cross-validation type
% ---------------------------------------------------------------------
cv_type        = cfg_choice;
cv_type.tag    = 'cv_type';
cv_type.name   = 'Cross-validation type';
cv_type.values = {cv_loso, cv_lkso, cv_losgo,cv_lksgo, cv_lobo,...
    cv_lkbo, cv_loro, cv_custom};
cv_type.val    = {cv_loso};
cv_type.help   = {'Choose the type of cross-validation to be used'};

% ---------------------------------------------------------------------
% machine_func Filename(s) of data
% ---------------------------------------------------------------------
machine_func        = cfg_files;
machine_func.tag    = 'machine_func';
machine_func.name   = 'Function';
machine_func.ufilter = '^*.m';
machine_func.num    = [1 1];
machine_func.help   = {'Choose a function that will perform prediction.'};

% ---------------------------------------------------------------------
% machine_args Regression Targets
% ---------------------------------------------------------------------
machine_args         = cfg_entry;
machine_args.tag     = 'machine_args';
machine_args.name    = 'Arguments';
machine_args.help    = {['Arguments for prediction machine.']};
machine_args.strtype = 's';
machine_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% custom_machine Regression group
% ---------------------------------------------------------------------
custom_machine         = cfg_branch;
custom_machine.tag     = 'custom_machine';
custom_machine.name    = 'Custom machine';
custom_machine.help    = {'Choose another prediction machine'};
custom_machine.val     = {machine_func, machine_args};

% ---------------------------------------------------------------------
% svm_opt SVM : flag whether to optimize C
% ---------------------------------------------------------------------
svm_opt         = cfg_menu;
svm_opt.tag     = 'svm_opt';
svm_opt.name    = 'Optimize hyper-parameter';
svm_opt.help    = {['Whether to optimize C, the SVM hyper-parameter, or not. '...
    'If Yes, than provide a range of possible values for C, in the form '...
    'min:step:max. Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100. ' ...
    'If not, a default value will be used (C=1).']};
svm_opt.labels  = {
    'No'
    'Yes'
    }';
svm_opt.values  = {0 1};
svm_opt.val     = {0};

% ---------------------------------------------------------------------
% svm_args Regression Targets
% ---------------------------------------------------------------------
svm_args         = cfg_entry;
svm_args.tag     = 'svm_args';
svm_args.name    = 'Soft-margin hyper-parameter';
svm_args.help    = {['Value(s) for prt_machine_svm_bin: soft-margin C. ',...
    'Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100.']};
svm_args.strtype = 'e';
svm_args.val     = {def.model.svmargs};
svm_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% cv_type Cross-validation type
% ---------------------------------------------------------------------
svm_cv_type_nested        = cfg_choice;
svm_cv_type_nested.tag    = 'cv_type_nested';
svm_cv_type_nested.name   = 'Cross-validation type for hyper-parameter optimization';
svm_cv_type_nested.values = {cv_loso,cv_lkso, cv_losgo,cv_lksgo, cv_lobo,...
    cv_lkbo, cv_loro};
svm_cv_type_nested.val    = {cv_loso};
svm_cv_type_nested.help   = {'Choose the type of cross-validation to be used'};

% ---------------------------------------------------------------------
% svm group
% ---------------------------------------------------------------------
svm         = cfg_branch;
svm.tag     = 'svm';
svm.name    = 'SVM Classification';
svm.help    = {'Binary support vector machine.'};
svm.val     = {svm_opt, svm_args, svm_cv_type_nested};

% ---------------------------------------------------------------------
% gpc_args GPC arguments
% ---------------------------------------------------------------------
gpc_args         = cfg_entry;
gpc_args.tag     = 'gpc_args';
gpc_args.name    = 'Arguments';
gpc_args.help    = {['Arguments for prt_machine_gpml']};
gpc_args.strtype = 's';
gpc_args.val     = {def.model.gpcargs};
gpc_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% gpc GPC
% ---------------------------------------------------------------------
gpc         = cfg_branch;
gpc.tag     = 'gpc';
gpc.name    = 'Gaussian Process Classification';
gpc.help    = {'Gaussian Process Classification'};
gpc.val     = {gpc_args};

% ---------------------------------------------------------------------
% gpclap_args GPC arguments
% ---------------------------------------------------------------------
gpclap_args         = cfg_entry;
gpclap_args.tag     = 'gpclap_args';
gpclap_args.name    = 'Arguments';
gpclap_args.help    = {['Arguments for prt_machine_gpclap']};
gpclap_args.strtype = 's';
gpclap_args.val     = {def.model.gpclapargs};
gpclap_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% gpclap GPC
% ---------------------------------------------------------------------
gpclap         = cfg_branch;
gpclap.tag     = 'gpclap';
gpclap.name    = 'Multiclass GPC';
gpclap.help    = {'Multiclass GPC'};
gpclap.val     = {gpclap_args};

% ---------------------------------------------------------------------
% sMKL_cla_opt L1-MKL : flag whether to optimize C
% ---------------------------------------------------------------------
sMKL_cla_opt         = cfg_menu;
sMKL_cla_opt.tag     = 'sMKL_cla_opt';
sMKL_cla_opt.name    = 'Optimize hyper-parameter';
sMKL_cla_opt.help    = {['Whether to optimize C, the SVM hyper-parameter, or not. '...
    'If Yes, than provide a range of possible values for C, in the form '...
    'min:step:max. Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100. ' ...
    'If not, a default value will be used (C=1).']};
sMKL_cla_opt.labels  = {
    'No'
    'Yes'
    }';
sMKL_cla_opt.values  = {0 1};
sMKL_cla_opt.val     = {0};

% ---------------------------------------------------------------------
% sMKL_cla_args L1-MKL arguments
% ---------------------------------------------------------------------
sMKL_cla_args         = cfg_entry;
sMKL_cla_args.tag     = 'sMKL_cla_args';
sMKL_cla_args.name    = 'Arguments';
sMKL_cla_args.help    = {['Arguments for prt_machine_sMKL_cla (same as for SVM)',...
    'Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100.']};
sMKL_cla_args.strtype = 'e';
sMKL_cla_args.val     = {def.model.l1MKLargs};
sMKL_cla_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% cv_type Cross-validation type
% ---------------------------------------------------------------------
sMKL_cla_cv_type_nested        = cfg_choice;
sMKL_cla_cv_type_nested.tag    = 'cv_type_nested';
sMKL_cla_cv_type_nested.name   = 'Cross-validation type for hyper-parameter optimization';
sMKL_cla_cv_type_nested.values = {cv_loso,cv_lkso, cv_losgo,cv_lksgo, cv_lobo,...
    cv_lkbo, cv_loro,cv_custom};
sMKL_cla_cv_type_nested.val    = {cv_loso};
sMKL_cla_cv_type_nested.help   = {'Choose the type of cross-validation to be used'};

% ---------------------------------------------------------------------
% sMKL_cla simple (L1) MKL
% ---------------------------------------------------------------------
sMKL_cla         = cfg_branch;
sMKL_cla.tag     = 'sMKL_cla';
sMKL_cla.name    = 'L1 Multi-Kernel Learning';
sMKL_cla.help    = {'Multi-Kernel Learning. Choose only if multiple kernels ' ...
    'were built during the feature set construction (either multiple modalities or ROIs). ' ...
    'It is strongly advised to "normalize" the kernels (in "operations").'};
sMKL_cla.val     = {sMKL_cla_opt, sMKL_cla_args, sMKL_cla_cv_type_nested};

% ---------------------------------------------------------------------
% gpr_args GPR arguments
% ---------------------------------------------------------------------
gpr_args         = cfg_entry;
gpr_args.tag     = 'gpr_args';
gpr_args.name    = 'Arguments';
gpr_args.help    = {['Arguments for prt_machine_gpr']};
gpr_args.strtype = 's';
gpr_args.val     = {def.model.gprargs};
gpr_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% gpr GPR
% ---------------------------------------------------------------------
gpr         = cfg_branch;
gpr.tag     = 'gpr';
gpr.name    = 'Gaussian Process Regression';
gpr.help    = {'Gaussian Process Regression'};
gpr.val     = {gpr_args};

% ---------------------------------------------------------------------
% sMKL_reg_opt L1-MKL : flag whether to optimize the hyperparameter
% ---------------------------------------------------------------------
sMKL_reg_opt         = cfg_menu;
sMKL_reg_opt.tag     = 'sMKL_reg_opt';
sMKL_reg_opt.name    = 'Optimize hyper-parameter';
sMKL_reg_opt.help    = {['Whether to optimize C, the MKL hyper-parameter, or not. '...
    'If Yes, than provide a range of possible values for C, in the form '...
    'min:step:max. Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100. ' ...
    'If not, a default value will be used (C=1).']};
sMKL_reg_opt.labels  = {
    'No'
    'Yes'
    }';
sMKL_reg_opt.values  = {0 1};
sMKL_reg_opt.val     = {0};

% ---------------------------------------------------------------------
% sMKL_reg_args sMKL_reg arguments
% ---------------------------------------------------------------------
sMKL_reg_args         = cfg_entry;
sMKL_reg_args.tag     = 'sMKL_reg_args';
sMKL_reg_args.name    = 'Arguments';
sMKL_reg_args.help    = {['Arguments for prt_machine_sMKL_reg']};
sMKL_reg_args.strtype = 'e';
sMKL_reg_args.val     = {def.model.l1MKLargs};
sMKL_reg_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% cv_type Cross-validation type
% ---------------------------------------------------------------------
sMKL_reg_cv_type_nested        = cfg_choice;
sMKL_reg_cv_type_nested.tag    = 'cv_type_nested';
sMKL_reg_cv_type_nested.name   = 'Cross-validation type for hyper-parameter optimization';
sMKL_reg_cv_type_nested.values = {cv_loso,cv_lkso, cv_losgo,cv_lksgo, cv_lobo,...
    cv_lkbo, cv_loro,cv_custom};
sMKL_reg_cv_type_nested.val    = {cv_loso};
sMKL_reg_cv_type_nested.help   = {'Choose the type of cross-validation to be used'};

% ---------------------------------------------------------------------
% sMKL_reg sMKL regression
% ---------------------------------------------------------------------
sMKL_reg         = cfg_branch;
sMKL_reg.tag     = 'sMKL_reg';
sMKL_reg.name    = 'Multi-Kernel Regression';
sMKL_reg.help    = {'Multi-Kernel Regression'};
sMKL_reg.val     = {sMKL_reg_opt, sMKL_reg_args, sMKL_reg_cv_type_nested};

% ---------------------------------------------------------------------
% krr_opt SVM : flag whether to optimize C
% ---------------------------------------------------------------------
krr_opt         = cfg_menu;
krr_opt.tag     = 'krr_opt';
krr_opt.name    = 'Optimize hyper-parameter';
krr_opt.help    = {['Whether to optimize K, the KRR hyper-parameter, or not. '...
    'If Yes, than provide a range of possible values for K, in the form '...
    'min:step:max. Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100. ' ...
    'If not, a default value will be used.']};
krr_opt.labels  = {
    'No'
    'Yes'
    }';
krr_opt.values  = {0 1};
krr_opt.val     = {0};

% ---------------------------------------------------------------------
% krr_args Regression Targets
% ---------------------------------------------------------------------
krr_args         = cfg_entry;
krr_args.tag     = 'krr_args';
krr_args.name    = 'Regularization';
krr_args.help    = {['Regularization for prt_machine_krr. ',...
    'Examples: 10.^[-2:5] or 1:100:1000 or 0.01 0.1 1 10 100.']};
krr_args.strtype = 'e';
krr_args.val     = {1};
krr_args.num     = [1 Inf];

% ---------------------------------------------------------------------
% cv_type Cross-validation type
% ---------------------------------------------------------------------
krr_cv_type_nested        = cfg_choice;
krr_cv_type_nested.tag    = 'cv_type_nested';
krr_cv_type_nested.name   = 'Cross-validation type for hyper-parameter optimization';
krr_cv_type_nested.values = {cv_loso,cv_lkso, cv_losgo,cv_lksgo, cv_lobo,...
    cv_lkbo, cv_loro,cv_custom};
krr_cv_type_nested.val    = {cv_loso};
krr_cv_type_nested.help   = {'Choose the type of cross-validation to be used'};

% ---------------------------------------------------------------------
% KRR group
% ---------------------------------------------------------------------
krr         = cfg_branch;
krr.tag     = 'krr';
krr.name    = 'Kernel Ridge Regression';
krr.help    = {'Kernel Ridge Regression.'};
krr.val     = {krr_opt,krr_args,krr_cv_type_nested};

% ---------------------------------------------------------------------
% RVR group
% ---------------------------------------------------------------------
rvr         = cfg_branch;
rvr.tag     = 'rvr';
rvr.name    = 'Relevance Vector Regression';
rvr.help    = {'Relevance Vector Regression. Tipping, Michael E.; Smola, Alex (2001).' ...
    '"Sparse Bayesian Learning and the Relevance Vector Machine". Journal of Machine Learning Research 1: 211?244.'};

% ---------------------------------------------------------------------
% rt_args Arguments to RT
% ---------------------------------------------------------------------
rt_args         = cfg_entry;
rt_args.tag     = 'rt_args';
rt_args.name    = 'Ntrees';
rt_args.help    = {['Number of trees in the forest.']};
rt_args.strtype = 'e';
rt_args.val     = {601};
rt_args.num     = [1 1];

% ---------------------------------------------------------------------
% RT group
% ---------------------------------------------------------------------
rt         = cfg_branch;
rt.tag     = 'rt';
rt.name    = 'Random Forest';
rt.help    = {'Random Forest. Breiman, Leo (2001)."Random Forests". ' ...
               'Machine Learning 45:5-32. This is a wrapper around ' ...
               'Peter Geurt''s implementation in his Regression Tree ' ...
               ' package.' };
rt.val     = {rt_args};

% % ---------------------------------------------------------------------
% % machine Select Features
% % ---------------------------------------------------------------------
% machine        = cfg_choice;
% machine.tag    = 'machine';
% machine.name   = 'Machine';
% machine.values = {svm,gpc,krr,rvr,rt,custom_machine};
% machine.val    =  {svm};
% machine.help   = {...
%     ['Choose a prediction machine for this model']};

% ---------------------------------------------------------------------
% machine_cl Select Machine
% ---------------------------------------------------------------------
machine_cl       = cfg_choice;
machine_cl.tag    = 'machine_cl';
machine_cl.name   = 'Machine';
machine_cl.values = {svm, gpc, gpclap, sMKL_cla, custom_machine}; 
machine_cl.val    =  {svm};
machine_cl.help   = {...
    ['Choose a prediction machine for this model']};
% Random Trees out since only kernel methods available for the moment
% machine_cl.values = {svm,gpc,gpclap,rt,sMKL_cla,custom_machine};

% ---------------------------------------------------------------------
% machine_rg Select Machine
% ---------------------------------------------------------------------
machine_rg       = cfg_choice;
machine_rg.tag    = 'machine_rg';
machine_rg.name   = 'Machine';
machine_rg.values = {krr,rvr,gpr,sMKL_reg,custom_machine};
machine_rg.val    =  {krr};
machine_rg.help   = {...
    ['Choose a prediction machine for this model']};

% ---------------------------------------------------------------------
% regression Regression
% ---------------------------------------------------------------------
reggroups         = cfg_repeat;
reggroups.tag     = 'reggroups';
reggroups.name    = 'Groups';
reggroups.help    = {['Add one group to this regression model. Click ''new'' '...
                    'or ''repeat'' to add another group.']};
reggroups.num     = [1 Inf];
reggroups.values  = {reg_group};

% ---------------------------------------------------------------------
% regression Regression
% ---------------------------------------------------------------------
regression         = cfg_branch;
regression.tag     = 'regression';
regression.name    = 'Regression';
regression.help    = {'Add group data and machine for regression.'};
regression.val     = {reggroups, machine_rg};

% ---------------------------------------------------------------------
% classes Classes
% ---------------------------------------------------------------------
classes         = cfg_repeat;
classes.tag     = 'classes';
classes.name    = 'Classes';
classes.help    = {['Specify which elements belong to this class. Click ''new'' '...
                           'or ''repeat'' to add another class.']};
classes.num     = [1 Inf];
classes.values  = {class};

% ---------------------------------------------------------------------
% classification Classification
% ---------------------------------------------------------------------
classification         = cfg_branch;
classification.tag     = 'classification';
classification.name    = 'Classification';
classification.help    = {'Specify classes and machine for classification.'};
classification.val     = {classes, machine_cl};

% ---------------------------------------------------------------------
% model_type Model type
% ---------------------------------------------------------------------
model_type        = cfg_choice;
model_type.tag    = 'model_type';
model_type.name   = 'Model Type ';
model_type.values = {classification, regression};
model_type.help   = {'Select which kind of predictive model is to be used.'};

% ---------------------------------------------------------------------
% include_allscans Include unused scans
% ---------------------------------------------------------------------
include_allscans         = cfg_menu;
include_allscans.tag     = 'include_allscans';
include_allscans.name    = 'Include all scans';
include_allscans.labels  = {
    'Yes'
    'No'
}';
include_allscans.values  = {1 0};
include_allscans.val     = {0};
include_allscans.help    = {[...
    'This option can be used to pass all the scans for each subject to ',...
    'the learning machine, regardless of whether they are directly ',...
    'involved in the classification or regression problem. For example, ',...
    'this can be used to estimate a GLM from the whole timeseries ',...
    'for each subject prior to prediction. This would allow the resulting ',...
    'regression coefficient images to be used as samples.']};

% ---------------------------------------------------------------------
% no_op All scans
% ---------------------------------------------------------------------
no_op         = cfg_const;
no_op.tag     = 'no_op';
no_op.name    = 'No operations';
no_op.val     = {1};
no_op.help    = {['No design specified. This option can be used '...
    'for modalities (e.g. structural scans) that do not '...
    'have an experimental design or for an fMRI design',...
    'where you want to include all scans in the timeseries']};

% ---------------------------------------------------------------------
% data_op Operation
% ---------------------------------------------------------------------
data_op         = cfg_menu;
data_op.tag     = 'data_op';
data_op.name    = 'Operation';
data_op.help    = {'Select an operation to apply.'};
% data_op.labels  = {
%     'Done'
%     'Sample averaging (within block)'
%     'Sample averaging (within subject/condition)'
%     'Mean centre features using training data'
%     'Normalize samples'
%     'Regress out covariates (subjects only)'
% }';
% data_op.values  = {0 1 2 3 4 5};
data_op.labels  = {
    'Done'
    'Sample averaging (within block)'
    'Sample averaging (within subject/condition)'
    'Mean centre features using training data'
    'Normalize samples'
    }';
data_op.values  = {0 1 2 3 4};
data_op.val     = {0};

% ---------------------------------------------------------------------
% data_op Operation
% ---------------------------------------------------------------------
data_op_mc         = cfg_menu;
data_op_mc.tag     = 'data_op_mc';
data_op_mc.name    = 'Mean centre features';
data_op_mc.help    = {'Select an operation to apply.'};
data_op_mc.labels  = {
    'Yes'
    'No'
}';
data_op_mc.values  = {1 0};
data_op_mc.val     = {1};

% ---------------------------------------------------------------------
% other_ops Other Operations
% ---------------------------------------------------------------------
other_ops         = cfg_repeat;
other_ops.tag     = 'other_ops';
other_ops.name    = 'Select Operations';
other_ops.help    = {...
    ['Add zero or more operations to be applied to the data before the ',...
     'prediction machine is called. These are executed within the ',...
     'cross-validation loop (i.e. they respect training/test independence) ',...
     'and will be executed in the order specified. ']};
other_ops.num     = [1 Inf];
other_ops.values  = {data_op};

% ---------------------------------------------------------------------
% use_other_ops Use other operations
% ---------------------------------------------------------------------
use_other_ops        = cfg_choice;
use_other_ops.tag    = 'use_other_ops';
use_other_ops.name   = 'Other Operations';
use_other_ops.values = {no_op, other_ops };
use_other_ops.val    = {no_op};
use_other_ops.help   = {'Include other operations?'};

% ---------------------------------------------------------------------
% sel_ops Class
% ---------------------------------------------------------------------
sel_ops         = cfg_branch;
sel_ops.tag     = 'sel_ops';
sel_ops.name    = 'Data operations';
sel_ops.help    = {...
    ['Specify operations to apply']};
sel_ops.val     = {data_op_mc use_other_ops};

% ---------------------------------------------------------------------
% data_ops Select Features
% ---------------------------------------------------------------------
data_ops        = cfg_choice;
data_ops.tag    = 'data_ops';
data_ops.name   = 'Data Operations';
data_ops.values = {data_op_mc, sel_ops};
data_ops.val    =  {no_op};
data_ops.help   = {...
    ['This branch controls operations that can be applied to the data ',...
     'before the data is passed to the classifier. Add zero or more ',...
     'operations to be applied. These will be executed in the order ',...
     'specified. ']}; 
 
% ---------------------------------------------------------------------
% model Model
% ---------------------------------------------------------------------
model        = cfg_exbranch;
model.tag    = 'model';
model.name   = 'Specify model';
model.val    = {infile, ...
                model_name, ...
                use_kernel, ...
                fsets, ...
                model_type, ...
                cv_type,...
                include_allscans,...
                sel_ops};
model.help   = {'Construct model according to design specified'};
model.prog   = @prt_run_model;
model.vout   = @vout_data;

%------------------------------------------------------------------------
%% Output function
%------------------------------------------------------------------------
function cdep = vout_data(job)
% Specifies the output from this modules, i.e. the filename of the mat file

cdep(1)            = cfg_dep;
cdep(1).sname      = 'PRT.mat file';
cdep(1).src_output = substruct('.','files');
cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
cdep(2)            = cfg_dep;
cdep(2).sname      = 'Model name';
cdep(2).src_output = substruct('.','mname');
cdep(2).tgt_spec   = cfg_findspec({{'strtype','s'}});
%------------------------------------------------------------------------

