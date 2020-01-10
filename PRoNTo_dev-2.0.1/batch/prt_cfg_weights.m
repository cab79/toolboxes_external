function weights = prt_cfg_weights
% Preprocessing of the data.
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M.J.Rosa
% $Id$

% ---------------------------------------------------------------------
% filename Filename(s) of data
% ---------------------------------------------------------------------
infile        = cfg_files;
infile.tag    = 'infile';
infile.name   = 'Load PRT.mat';
infile.filter = 'mat';
infile.num    = [1 1];
infile.help   = {['Select PRT.mat (file containing data/design/models ',...
                 'structure).']};

% ---------------------------------------------------------------------
% model_name Feature set name
% ---------------------------------------------------------------------
model_name         = cfg_entry;
model_name.tag     = 'model_name';
model_name.name    = 'Model name';
model_name.help    = {'Name of a model. Must correspond with one ' ...
					  'specified in ''Run model'' and ''Specify model'''...
                      'batch modules.' };
model_name.strtype = 's';
model_name.num     = [1 Inf];

% ---------------------------------------------------------------------
% img_name Feature set name
% ---------------------------------------------------------------------
img_name         = cfg_entry;
img_name.tag     = 'img_name';
img_name.name    = 'Image name (optional)';
img_name.help    = {['Name of the file with weights (optional). If left empty ',...
                    ' an automatic name will be generated.']};
img_name.strtype = 's';
img_name.num     = [0 Inf];
img_name.val     = {''};

% ---------------------------------------------------------------------
% flag_cwi Build the weight images for each permutation (optional)
% ---------------------------------------------------------------------
flag_cwi         = cfg_menu;
flag_cwi.tag     = 'flag_cwi';
flag_cwi.name    = 'Build weight images for permutations';
flag_cwi.help    = {['Set to Yes to compute the weight images obtained ' ...
    'from each permutation. This is to further assess the significance ' ...
    'of the ranking distance between two models.']};
flag_cwi.labels  = {
               'Yes'
               'No'
}';
flag_cwi.values  = {1 0};
flag_cwi.val     = {0};

% ---------------------------------------------------------------------
% atl_name Filename for the atlas to compute the weights per ROI
% ---------------------------------------------------------------------
atl_name         = cfg_files;
atl_name.tag     = 'atl_name';
atl_name.name    = 'Load Atlas';
atl_name.ufilter = '.*';
atl_name.filter  = 'image';
atl_name.num     = [0 1];
atl_name.val     = {{''}};
atl_name.def     = @(val)prt_get_defaults('fs.atlasroi', val{:});
atl_name.help    = {'Select atlas file to build weights per ROI.'};

% ---------------------------------------------------------------------
% no_atl No weight per ROI
% ---------------------------------------------------------------------
no_atl         = cfg_const;
no_atl.tag     = 'no_atl';
no_atl.name    = 'No weight per ROI';
no_atl.val     = {0};
no_atl.help    = {'Not computing weight per ROI.'};

% ---------------------------------------------------------------------
% build_wpr Build the weight images per ROI
% ---------------------------------------------------------------------
build_wpr         = cfg_choice;
build_wpr.tag     = 'build_wpr';
build_wpr.name    = 'Build weight images per ROI';
build_wpr.help    = {['Set to Yes to compute the weight images per ROI ' ...
    'You need then to select the atlas image.']};
% build_wpr.labels  = {
%                'Yes'
%                'No'
% }';
build_wpr.values  = {no_atl atl_name};
build_wpr.val     = {no_atl};

% ---------------------------------------------------------------------
% flag_wroi Build the weight images for each ROI (optional)
% ---------------------------------------------------------------------
flag_wroi         = cfg_menu;
flag_wroi.tag     = 'flag_cwi';
flag_wroi.name    = 'Build weight images per ROI';
flag_wroi.help    = {['Set to Yes to compute the weight images obtained ' ...
    'from each ROI.']};
flag_wroi.labels  = {
               'Yes'
               'No'
}';
flag_wroi.values  = {1 0};
flag_wroi.val     = {0};

% ---------------------------------------------------------------------
% cv_model Preprocessing
% ---------------------------------------------------------------------
weights        = cfg_exbranch;
weights.tag    = 'weights';
weights.name   = 'Compute weights';
weights.val    = {infile model_name img_name build_wpr flag_cwi};
weights.help   = {[
    'Compute weights. This module computes the linear weights of a classifier ',...
    'and saves them as a 4D image. 3 dimensions correspond to the image dimensions specified in ',...
    'the second-level mask, while the extra dimension corresponds to the number of folds. ',...
    'There is one 3D weights image per fold.']};
weights.prog   = @prt_run_weights;
weights.vout   = @vout_data;

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
cdep(2).sname      = 'Weight image file';
cdep(2).src_output = substruct('.','files');
cdep(2).tgt_spec   = cfg_findspec({{'filter','img','strtype','e'}});
%------------------------------------------------------------------------
