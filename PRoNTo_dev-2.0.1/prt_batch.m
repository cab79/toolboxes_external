function prt_batch
% Pattern Recognition for Neuroimaging Toolbox, PRoNTo.
%
% This function prepares and launches the batch system.
% This builds the whole tree for the various tools and their GUI at the
% first call to this script.
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

persistent batch_initialize
global PRT_INIT

if isempty(PRT_INIT) || ~PRT_INIT
    prt('startup','nogui');
end

if isempty(batch_initialize) || ~batch_initialize
    % PRoNTo config tree
    prt_gui = prt_cfg_batch;
    % Adding PRoNTo config tree to the SPM tools
    cfg_util('addapp', prt_gui)
    % No need to do it again for this session
    batch_initialize = 1;
end

% Launching the batch system
cfg_ui

return
