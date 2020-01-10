function varargout = prt_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
%
% FORMAT defaults = prt_get_defaults
% Return the global "defaults" variable defined in prt_defaults.m.
%
% FORMAT defval = prt_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr". 
% Currently, this is a '.' subscript reference into the global  
% "prt_def" variable defined in prt_defaults.m.
%
% FORMAT prt_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of PRoNTo. To make
% persistent changes, edit prt_defaults.m.
%
% The structure and content of this file are largely inspired by SPM &
% Matlabbatch.
% http://www.fil.ion.ucl.ac.uk/spm
% http://sourceforge.net/projects/matlabbatch/
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Originally written by Volkmar Glauche
% Then modified for use with the PRoNTo toolbox by Christophe Phillips
% $Id$

global prt_def;
if isempty(prt_def)
    prt_defaults;
end

if nargin == 0
    varargout{1} = prt_def;
    return
end

% construct subscript reference struct from dot delimited tag string
tags = textscan(defstr,'%s', 'delimiter','.');
subs = struct('type','.','subs',tags{1}');

if nargin == 1
    varargout{1} = subsref(prt_def, subs);
else
    prt_def = subsasgn(prt_def, subs, varargin{1});
end
