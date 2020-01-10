function varargout = prt(varargin)
% Pattern Recognition for Neuroimaging Toolbox, PRoNTo.
%
% This function initializes things for PRoNTo and provides some low level
% functionalities
%
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

% TODO:
% - fix which subdirectories from all the machines are necessary, and only
%   add these to Matlab path.

%-Format arguments
%-----------------------------------------------------------------------
global PRT_INIT
if nargin == 0,
    Action = 'StartUp';
else
    Action = varargin{1};
end

switch lower(Action)
    %==================================================================
    case 'startup'                                    % Startup the toolbox
        %==================================================================
        
        % Welcome message
        prt('ASCIIwelcome');
        
        % add appropriate paths, if necessary
        %   - batch dir
        if ~exist('prt_cfg_batch','file')
            addpath(fullfile(prt('Dir'),'batch'));
        end
        %   - machines
        if ~exist('prt_machine','file')
            pth_machines = fullfile(prt('Dir'),'machines');
            addpath(pth_machines);
            % add each machine's sub-directory
            % and ALL its subdirectories recursively
            ls_machinedir = list_subdir(pth_machines);
            for ii=1:numel(ls_machinedir)
                gpath_ii = genpath(fullfile(pth_machines,ls_machinedir{ii}));
                gpath_ii = clean_gpath(gpath_ii);
                addpath(gpath_ii)
            end
        end
        
        % utils - dirty check for the moment
        if ~exist('prt_checkAlphaNumUnder','file')
            addpath(fullfile(prt('Dir'),'utils'));
        end
        
        % check installation of machines and that of SPM8/12
        ok = check_installation;
        if ~ok
            beep
            fprintf('INSTALLATION PROBLEM!');
            return
        end
        
        % Add SPM's directories: matlabbatch
        if ~exist('cfg_util','file')
            addpath(fullfile(spm('Dir'),'matlabbatch'));
        end
        
        % intialize the matlabbatch system
        cfg_get_defaults('cfg_util.genscript_run', @genscript_run);
        cfg_util('initcfg');
        clear prt_batch;
        
        % set path to PRoNTo and SPM dir into 'file select'
        spm_select('prevdirs',[spm('Dir') filesep]);
        spm_select('prevdirs',[prt('Dir') filesep]);
        
        % launch the main GUI, if needed
        if nargin<2 || ~strcmp(varargin{2},'nogui')
            prt_ui_main;
        end
        
        % print present working directory
        fprintf('PRoNTo present working directory:\n\t%s\n',pwd)
        
        % Init flag true
        PRT_INIT = true;
        
        %==================================================================
    case 'asciiwelcome'                       %-ASCII PRoNTo banner welcome
        %==================================================================
        disp( '                                                             ');
        disp( '     ____  ____        _   ________              ___    ____ ');
        disp( '    / __ \/ __ \____  / | / /_  __/___     _   _|__ \  / __ \');
        disp( '   / /_/ / /_/ / __ \/  |/ / / / / __ \   | | / /_/ / / / / /');
        disp( '  / ____/ _, _/ /_/ / /|  / / / / /_/ /   | |/ / __/_/ /_/ / ');
        disp( ' /_/   /_/ |_|\____/_/ |_/ /_/  \____/    |___/____(_)____/  ');
        disp( '                                                             ');
        disp( '      PRoNTo v2.0 - http://www.mlnl.cs.ucl.ac.uk/pronto      ');
        fprintf('\n');
        
        %==================================================================
    case 'dir'                          %-Identify specific (PRT) directory
        %==================================================================
        % prt('Dir',Mfile)
        %------------------------------------------------------------------
        if nargin<2,
            Mfile = 'prt';
        else
            Mfile = varargin{2};
        end
        PRTdir = which(Mfile);
        
        if isempty(PRTdir)             %-Not found or full pathname given
            if exist(Mfile,'file')==2  %-Full pathname
                PRTdir = Mfile;
            else
                error(['Can''t find ',Mfile,' on MATLABPATH']);
            end
        end
        PRTdir    = fileparts(PRTdir);
        varargout = {PRTdir};
 
        %==================================================================
    case 'ver'                                             %-PRoNTo version
        %==================================================================
        % [ver, rel] = prt('Ver',ReDo)
        %------------------------------------------------------------------
        % NOTE:
        % 1/
        % This bit of code is largely inspired/copied from SPM8!
        % See http://www.fil.ion.ucl.ac.uk/spm for details.
        %
        % 2/
        % NOT usable any more to get the version of an individual file 
        % since the  switch from SVN to GitHub!
        % Git does not update an Id tag inside the files, hence no option 
        % to check the version automatically...
       
        if nargin ~= 2,
            ReDo = false;
        else
            ReDo = logical(varargin{3});
        end
        
        v = get_version(ReDo);
        varargout = {v.Release v.Version};
        
        %==================================================================
    otherwise                                       %-Unknown action string
        %==================================================================
        error('Unknown action string');
        
end

return

%=======================================================================
%% SUBFUNCTIONS
%=======================================================================

%=======================================================================
function ok = check_installation
%=======================================================================
% Function to check installation state of machines and SPM

ok = true;

% Check SPM installation
if exist('spm.m','file')
    [SPMver, SPMrel] = spm('Ver');
    if (~(strcmpi(SPMver,'spm8') && str2double(SPMrel)>8.5)) && ...
            isempty(regexpi(SPMver,'spm12'))
        beep
        fprintf('\nERROR:\n')
        fprintf('\tThe *latest* version of SPM8 or SPM12 should be installed on your computer,\n')
        fprintf('\tand be available on MATLABPATH!\n\n')
        ok = false;
    end
else
    beep
    fprintf('\nERROR:\n')
    fprintf('\tThe *latest* version of SPM8 or SPM12 should be installed on your computer,\n')
    fprintf('\tand be available on MATLABPATH!\n\n')
    ok = false;
end


% Check for compiled routines
%============================
% - svm
dumb = which('svmtrain');
if ~isempty(strfind(dumb,'libsvm'))          % svm found in libsvm folder, OK
    disp('SVM path: OK')
    flag=1;
elseif ~isempty(strfind(dumb,'biolearning')) % svm found in the Matlab toolbox
    flag=0;
else                                         % svm not found at all, so need to compile
    flag=2;
end

% in case PRoNTo was installed with all subfolders under the biostats
if ~flag
    pth_machines = fullfile(prt('Dir'),'machines');
    addpath(pth_machines);
    % add each machine's sub-directory
    % and ALL its subdirectories recursively
    ls_machinedir = list_subdir(pth_machines);
    for ii=1:numel(ls_machinedir)
        gpath_ii = genpath(fullfile(pth_machines,ls_machinedir{ii}));
        gpath_ii = clean_gpath(gpath_ii);
        addpath(gpath_ii)
    end
    dumb = which('svmtrain');
    if isempty(strfind(dumb,'libsvm'))
        flag=2; %s till not working, need to recompile
    elseif ~isempty(strfind(dumb,'biolearning'))
        flag = 2;
        disp('PRoNTo was found under the biostats toolbox, please correct path')
        disp('SVM path: OK')
    else
        flag =1;
    end
end

if flag ==2 % need to recompile for the OS
    pth_machines = fullfile(prt('Dir'),'machines');
    ls_machinedir = list_subdir(pth_machines);
    for i=1:length(ls_machinedir)
        if ~isempty(strfind(ls_machinedir{i},'libsvm'))
            pfn= fullfile(pth_machines,ls_machinedir{i});
            dirtorem=cd;
            cd(pfn)
            cd matlab
            make;
            cd(dirtorem)
        end
    end
    dumb = which('svmtrain');
    if isempty(strfind(dumb,'libsvm'))
        %could not recompile
        beep
        warning('PRoNTo:SVMcompilation', ...
        ['SVM path not recognized. Please check that: \n', ...
        '- PRoNTo''directory was added *without* all subfolders \n',...
        '- PRoNTo is above the biostats Matlab toolbox \n',...
        'Otherwise, the routines surely need to be re-compiled for your OS \n',...
        'Please look on the web or ask on the mailing list for assistance'])
    end
end

% - GP
dumb = which('solve_chol');
if ~isempty(dumb) && ~isempty(strfind(dumb,'.mex'))
    disp('GP path: OK')
else
    beep
    disp('GP not compiled: routines will work but be slower')
end

% NOTE:
% Tree-based methods not available in this version. 
% So no need to check for it!
% We still plan to have it in a future release.
% 
% % - RF
% dumb = which('rtenslearn_c');
% if ~isempty(dumb) || ~isempty(strfind(dumb,'.mex'))
%     disp('RF path: OK')
% else
%     beep
%     warning('PRoNTo:RFcompilation', ...
%         ['RF path not recognized. Please check that \n', ...
%         'PRoNTo was added with all subfolders \n',...
%         'Otherwise, the routines surely need to be re-compiled for your OS \n',...
%         'Please look on the web or ask on the mailing list for assistance'])
% end

return

%=======================================================================
function lsdir = list_subdir(pth_dir,rejd)
%=======================================================================
% function that returns the list of subdirectories of a directory,
% rejecting those beginning with some characters ('.', '@' and '_' by
% default)

if nargin<2
    rejd = '.@_';
end
if nargin<1
    pth_dir = pwd;
end

tmp = dir(pth_dir);
ld = find([tmp.isdir]); ld([1 2]) = [];
lsdir = {tmp(ld).name};
if ~isempty(rejd)
    for ii=1:numel(rejd)
        lrej = find(strncmp(rejd(ii),lsdir,1));
        if ~isempty(lrej)
            lsdir(lrej) = [];
        end
    end
end

return

%=======================================================================
function v = get_version(ReDo)                 %-Retrieve PRoNTo version
%=======================================================================
persistent PRoNTo_VER;
v = PRoNTo_VER;
if isempty(PRoNTo_VER) || (nargin > 0 && ReDo)
    v = struct('Name','','Version','','Release','','Date','');
    try
        vfile = fullfile(prt('Dir'),'prt_contents.m');
        fid = fopen(vfile,'rt');
        if fid == -1, error(str); end
        l1 = fgetl(fid); l2 = fgetl(fid);
        fclose(fid);
        l1 = strtrim(l1(2:end)); l2 = strtrim(l2(2:end));
        t  = textscan(l2,'%s','delimiter',' '); t = t{1};
        v.Name = l1; v.Date = t{4};
        v.Version = t{2}; v.Release = t{3}(2:end-1);
    catch %#ok<CTCH>
        error('PRoNTo:getversion', ...
            'Can''t obtain PRoNTo Revision information.');
    end
    PRoNTo_VER = v;
end

return

%=======================================================================
function gpath = clean_gpath(gpath,rejd)
%=======================================================================
% function that "cleans up" a list of pathes to subdirectories,
% i.e. it removes any path containing a set of strings.
% By default, it removes all the '.svn' pathes. Other strings can be passed
% as a cell array

if nargin<2
    rejd = {'.svn'};
end
if nargin<1
    return
end

if numel(rejd)>1
    % do it 1 by 1
    for ii=1:numel(rejd)
        gpath = clean_gpath(gpath,rejd{ii});
    end
else
    % deal with 1 string
    l_col = strfind(gpath,':');
    for ii=numel(l_col):-1:1
        if ii>1
            pth_bit = [l_col(ii-1)+1 l_col(ii)];
        else
            pth_bit = [1 l_col(ii)];
        end
        if ~isempty(strfind(gpath(pth_bit(1):pth_bit(2)),rejd{1}))
            % remove the bit
            gpath(pth_bit(1):pth_bit(2)) = [];
        end
    end
end

return
