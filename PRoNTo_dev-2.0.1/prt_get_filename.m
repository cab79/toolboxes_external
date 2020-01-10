function out = prt_get_filename(ids)

% out = prt_get_filename(ids)
% 
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand
% $Id$

% Obs: Names return to the simple form: e.g. 'g1_s1_m1_c1'

try
    group_prefix = ['g' int2str(ids(1)) '_'];
catch
    group_prefix = '';
end

try
    subj_prefix = ['s' int2str(ids(2)) '_'];
catch
    subj_prefix = '';
end

try
    mod_prefix = ['m' int2str(ids(3))];
catch
    mod_prefix = '';
end

try
    cond_prefix = ['_c' int2str(ids(4))];
catch
    cond_prefix = '';
end

out = ['PRT_' group_prefix subj_prefix mod_prefix cond_prefix];
