function out = prt_run_cv_model(varargin)
%
% PRONTO job execution function
%
% INPUT
%   job    - harvested job data structure (see matlabbatch help)
%
% OUTPUT
%   out    - filename of saved data structure (1 file per group, per 
%            subject, per modality, per condition 
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand
% $Id$

job   = varargin{1};

% Load PRT.mat
% -------------------------------------------------------------------------
fname = char(job.infile);
if exist('PRT','var')
    clear PRT
end
load(fname);

% -------------------------------------------------------------------------
% Input file
% -------------------------------------------------------------------------

in.fname      = job.infile;
in.model_name = job.model_name;
mid = prt_init_model(PRT, in);

% Special cross-validation for MCKR
if strcmp(PRT.model(mid).input.machine.function,'prt_machine_mckr')
    fname = prt_cv_mckr(PRT,in);
else
    fname = prt_cv_model(PRT, in);
end

% Permutation test, required.
load(fname) % reload updated PRT!
if isfield(job,'perm_test') % to ensure back compatibility with older batch
    if isfield(job.perm_test,'perm_t')
        if isfield(job.perm_test.perm_t,'flag_sw') %keep compatibility
            flag = job.perm_test.perm_t.flag_sw;
        else
            flag = 0;
        end
        prt_permutation(PRT, job.perm_test.perm_t.N_perm, mid, ...
            spm_str_manip(fname,'h'),flag);
    end
end

% -------------------------------------------------------------------------
% Function output
% -------------------------------------------------------------------------
out = []; %prevent warning of overwriting 'char' class
disp('Model execution complete.')
out.files{1} = in.fname{1};
disp('Done')

 return