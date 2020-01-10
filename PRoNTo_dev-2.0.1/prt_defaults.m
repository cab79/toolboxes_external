function prt_defaults
% Sets the defaults which are used by the Pattern Recognition for
% Neuroimaging Toolbox, aka. PRoNTo.
%
% FORMAT prt_defaults
%_______________________________________________________________________
%
% This file can be customised to any the site/person own setup.
% Individual users can make copies which can be stored on their own
% matlab path. Make sure your 'prt_defaults' is the first one found in the
% path. See matlab documentation for details on setting path.
%
% Care must be taken when modifying this file!
%
% The structure and content of this file are largely inspired by SPM:
% http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

%%
global prt_def

% Global defaults
% prt_loc = which('prt_batch');
% prt_def.global.install_dir = fileparts(prt_loc);
prt_def.global.install_dir = prt('dir');

% Default colors of the different elements
%-----------------------------------------------
% prt_def.color.bg1=[0 0.8 1];
% prt_def.color.bg2=[0.9,0.6,0.3];
% prt_def.color.fr=[1,0.5,0.7];
% prt_def.color.high=[0.2,0.2,0.8];

prt_def.color.bg1  = [0.83,0.83,0.83];
prt_def.color.bg2  = [0.88,0.88,0.88];
prt_def.color.fr   = [0.92,0.92,0.92];
prt_def.color.high = [0.8 0 0];
prt_def.color.black = [0 0 0];

% Parameters for the data and design
%-----------------------------------------------
prt_def.datad.hrfd = 0; % HRF delay in seconds
prt_def.datad.hrfw = 0; % HRF FWHM, used to compute the overlap between conditions

prt_def.prep.default_mask  = fullfile(prt('dir'),'masks', ...
                                'SPM_mask_noeyes.hdr');% default mask

% Preprocessing defaults
%------------------------------------------------
% memory limit for kernel/file arrays construction
prt_def.fs.mem_limit = 256*1024*1024;  % bytes of memory to use
prt_def.fs.writeraw  = 0;              % flag to write the data detrended (default) or raw (to set to 1).

% Default atlas for ROI defintion
prt_def.fs.atlasroi  = cellstr(fullfile(prt('dir'),'atlas', ...
    'aal_79x91x69.img')); 

% Design specification default
prt_def.dspec.use3 = [1 2];



% Specify model: Parameters of the different machines
%--------------------------------------------------
prt_def.model.svmargs     = 1;
prt_def.model.libsvmargs  = '-q -s 0 -t 4 -c ';
prt_def.model.gpcargs     = '-l erf -h';%-h 
prt_def.model.gpclapargs  = '-h'; %'-h';
prt_def.model.gprargs     = '-l gauss -h'; % -h
prt_def.model.krrargs     = 1;
prt_def.model.rtargs      = 601;
prt_def.model.l1MKLargs   = 1;
prt_def.model.l1MKLmaxitr = 250;

% Parralelization of the code
%--------------------------------------------------
prt_def.paral.allow     = false; % use (or not) 'parfor' loops
prt_def.paral.ncore     = 3;     % number of cores that can be used.

return
