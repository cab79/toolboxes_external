function prt = prt_cfg_batch
% Pattern Recognition for Neuroimaging Toolbox, PRoNTo.
% PRT configuration file
% This builds the whole tree for the various tools and their GUI.
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

% NOTE:
% Only de-comment the following line if
% - PRT is installed in SPM8/12 toolbox directory, for example
%   SPM12\toolbox\PRoNTo
% - and this directory isn't saved on Matlab path (which it shouldn't)
%
% if ~isdeployed, addpath(fullfile(spm('Dir'),'toolbox','PRoNTo')); end

% ---------------------------------------------------------------------
%% prt PRoNTo series of modules
% ---------------------------------------------------------------------
prt         = cfg_choice;
prt.tag     = 'prt';
prt.name    = 'PRoNTo';
prt.help    = {[...
    'This is the batch interface for PRoNTo, i.e. Pattern Recognition '...
    'for Neuroimaging Toolbox providing a GUI for the various tools.']
                  }';
prt.values  = {prt_cfg_design, ...
               prt_cfg_fs, ...
               prt_cfg_model, ...
               prt_cfg_cv_model,...
               prt_cfg_weights};
%------------------------------------------------------------------------

return