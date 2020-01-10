function out = prt_run_design(varargin)
%
% PRoNTo job execution function
% takes a harvested job data structure and rearranges data into PRT
% data structure, then saves PRT.mat file.
%
% INPUT
%   job    - harvested job data structure (see matlabbatch help)
%
% OUTPUT
%   out    - filename of saved data structure.
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory
%
% Written by M.J.Rosa
% $Id$


% Job variable
% -------------------------------------------------------------------------
job   = varargin{1};

% Back compatibility (at least trying to...)
% -------------------------------------------------------------------------
if isfield(job,'hrfover') % old setup
    job.fmri_des.hrfover = job.hrfover;
    job.fmri_des.hrfdel  = job.hrfdel;
end

% Directory
% -------------------------------------------------------------------------
fname   = 'PRT.mat';
fname   = fullfile(job.dir_name{1},fname);

% Number of group
% -------------------------------------------------------------------------
ngroup    = length(job.group);

% Masks
% -------------------------------------------------------------------------
nmasks     = length(job.mask);
for i = 1:nmasks
    mod_names{i}      = job.mask(i).mod_name;
    masks(i).mod_name = mod_names{i};
    masks(i).fname    = char(job.mask(i).fmask);
end

mod_names_uniq = unique(mod_names);

if nmasks ~= length(mod_names_uniq);
    out.files{1} = [];
    beep;
    sprintf('Names of mask modalities repeated! Please correct!')
    return
end

% Make PRT.mat
% -------------------------------------------------------------------------

% Data type
if isfield(job.group(1).select,'modality')
    % selection by "images" in a modality
    nmod_scans = length(job.group(1).select.modality);
    for g = 1:ngroup
        
        clear mod_names_mod
        nmod   = length(job.group(g).select.modality);
        nsub   = length(job.group(g).select.modality(1).subjects);
        
        % Check if the number of masks and conditions is the same
        if nmod ~= nmod_scans
            out.files{1} = [];
            beep;
            sprintf('Numbers of modalities in groups 1 and %d differ!',g)
            disp('Please correct!')
            return
        else
            if nmod ~= nmasks
                out.files{1} = [];
                beep;
                sprintf('Number of modalities in group %d different from number of masks!',g)
                disp('Please correct!')
                return               
            else
                % Modalities
                PRT.group(g).gr_name  = job.group(g).gr_name;
                % Subjects
                for s = 1:nsub,
                    subj_name = sprintf('S%d',s);
                    for m = 1:nmod,
                        modnm   = job.group(g).select.modality(m).mod_name;
                        ns      = length(job.group(g).select.modality(m).subjects);
                        if ~isempty(job.group(g).select.modality(m).rt_subj)
                                rt_subj = job.group(g).select.modality(m).rt_subj(:);
                                if length(rt_subj) ~= ns
                                    out.files{1} = [];
                                    beep
                                    sprintf('Number of regression targets must be the number of subjects/scans! ')
                                    disp('Please correct!')
                                    return
                                else
                                    PRT.group(g).subject(s).modality(m).rt_subj = rt_subj(s);
                                end
                        else
                            PRT.group(g).subject(s).modality(m).rt_subj = [];
                        end
                        if ~isempty(job.group(g).select.modality(m).covar{1})
                            try
                                load(char(job.group(g).select.modality(m).covar{1}));
                                if exist('R','var')
                                    if size(R,1)==ns
                                        PRT.group(g).subject(s).modality(m).covar  = R(s,:);
                                    else
                                        out.files{1} = [];
                                        beep
                                        sprintf('Number of covariates must be the number of subjects/scans! ')
                                        disp('Please correct!')
                                        return
                                    end
                                else
                                    out.files{1} = [];
                                    beep
                                    sprintf('Covariates file must contain ''R'' variable! ')
                                    disp('Please correct!')
                                    return
                                end
                            catch
                                beep
                                sprintf('Could not load %s file!',char(job.group(g).select.modality(m).covar{1}))
                                out.files{1} = [];
                                return
                            end
                        else
                            PRT.group(g).subject(s).modality(m).covar  = [];
                        end
                        mod_names_mod{m} = modnm;
                        if isempty(intersect(mod_names_uniq,modnm)),
                            out.files{1} = [];
                            beep
                            sprintf('Incorrect modality name %s for subject %d group %d! ',modnm,s,g)
                            disp('Please correct!')
                            return
                        end
                        if nsub ~= ns
                            out.files{1} = [];
                            beep
                            sprintf('Number of subjects in modality %d and 1 of group %d are different! ',m,g)
                            disp('Please correct!')
                            return
                        else
                            PRT.group(g).subject(s).subj_name            = subj_name;
                            PRT.group(g).subject(s).modality(m).mod_name = job.group(g).select.modality(m).mod_name;
                            PRT.group(g).subject(s).modality(m).design   = 0;
                            PRT.group(g).subject(s).modality(m).scans    = char(job.group(g).select.modality(m).subjects{s});
                        end
                    end
                    if nmod ~= length(unique(mod_names_mod));
                        out.files{1} = [];
                        beep;
                        sprintf('Names of modalities in group %d repeated! Please correct!',g)
                        return
                    end
                end
            end
        end
        PRT.group(g).hrfoverlap = job.fmri_des.hrfover;
        PRT.group(g).hrfdelay   = job.fmri_des.hrfdel;
    end
else
    % selection by subject
    for g = 1:ngroup,  
        
        nmod_subjs = length(job.group(1).select.subject{1});
        nsubj  = length(job.group(g).select.subject);
        nsubj1 = length(job.group(1).select.subject);
        
        if nsubj ~= nsubj1,
            disp('Warning: unbalanced groups.')
        end
        for j = 1:nsubj,
            clear mod_names_subj
            subj_name = sprintf('S%d',j);
            nmod = length(job.group(g).select.subject{j});
            % Check if the number of masks and conditions is the same
            if nmod ~= nmod_subjs
                out.files{1} = [];
                beep
                sprintf('Numbers of modalities in subjects 1 and %d from group %d differ!',j,g)
                disp('Please correct!')
                return
            else
                if nmod ~= nmasks
                    out.files{1} = [];
                    beep
                    sprintf('Number of modalities in group %d subject %d different from number of masks!',g,j)
                    disp('Please correct!')
                    return
                else
                    for k = 1:nmod,
                        modnm    = job.group(g).select.subject{j}(k).mod_name;
                        TR       = job.group(g).select.subject{j}(k).TR;
                        mod_names_subj{k} = modnm;
                        if isempty(intersect(mod_names_uniq,modnm)),
                            out.files{1} = [];
                            beep
                            sprintf('Incorrect modality name %s for subject %d group %d! ',modnm,j,g)
                            disp('Please correct!')
                            return
                        end
                        clear design
                        if isfield(job.group(g).select.subject{j}(k).design,'load_SPM')
                            % Load SPM.mat design
                            try
                                load(job.group(g).select.subject{j}(k).design.load_SPM{1});
                            catch
                                out.files{1} = [];
                                beep
                                disp('Could not load SPM.mat file!')
                                return
                            end
                            switch lower(SPM.xBF.UNITS)
                                case 'scans'
                                    unit   = 0;
                                case 'seconds'
                                    unit   = 1;
                                case 'secs'
                                    unit   = 1;
                            end        
                            nscans  = length(job.group(g).select.subject{j}(k).scans);
                            ncond   = length(SPM.Sess(1).U);
                            for c = 1:ncond
                                conds(c).cond_name = SPM.Sess(1).U(c).name{1}; %#ok<*AGROW>
                                conds(c).onsets    = SPM.Sess(1).U(c).ons;
                                conds(c).durations = SPM.Sess(1).U(c).dur;
                            end                        
                            checked_conds = prt_check_design(conds,TR,unit,job.fmri_des.hrfover,job.fmri_des.hrfdel);
                            design.conds  = checked_conds.conds;
                            design.stats  = checked_conds.stats;
                            design.TR     = TR;
                            design.unit   = unit;
                            maxcond       = max([design.conds(:).scans]);
                            if nscans>1 && nscans < maxcond
                                sprintf('Design of subject %d, group %d, modality %d, exceeds time series!',j,g,k)
                                disp('Corresponding events were discarded')
                                for l = 1:length(design.conds)
                                    ovser = find(design.conds(l).scans > nscans);
                                    inser = find(design.conds(l).scans <= nscans);
                                    design.conds(l).discardedscans = [design.conds(l).discardedscans, design.conds(l).scans(ovser)];
                                    design.conds(l).scans = design.conds(l).scans(inser);
                                    design.conds(l).blocks = design.conds(l).blocks(inser);
                                end
                            end
                        else
                            if isfield(job.group(g).select.subject{j}(k).design,'no_design')
                                % No design
                                design = 0;
                            else
                                % Manual design
                                nscans = length(job.group(g).select.subject{j}(k).scans);
                                unit   = job.group(g).select.subject{j}(k).design.new_design.unit;
                                % Create new design
                                if ~isempty(job.group(g).select.subject{j}(k).design.new_design.multi_conds{1})
                                    multi_fname = job.group(g).select.subject{j}(k).design.new_design.multi_conds{1};
                                    % Multiple conditions
                                    try
                                        load(multi_fname);
                                    catch
                                        beep
                                        sprintf('Could not load %s file!',multi_fname)
                                        out.files{1} = [];
                                        return
                                    end
                                    try
                                        multicond.names = names;
                                    catch
                                        beep
                                        disp('No "names" found in the .mat file, please select another file!')
                                        out.files{1} = [];
                                        return
                                    end
                                    try
                                        multicond.durations = durations;
                                    catch
                                        beep
                                        disp('No "durations" found in the .mat file, please select another file!')
                                        out.files{1} = [];
                                        return
                                    end
                                    try
                                        multicond.onsets = onsets;
                                    catch
                                        beep
                                        disp('No "onsets" found in the .mat file, please select another file!')
                                        out.files{1} = [];
                                        return
                                    end
                                    try
                                        multicond.rt_trial = rt_trial;
                                    catch
                                        multicond.rt_trial = cell(length(multicond.onsets),1);
                                    end
                                    for mc = 1:length(multicond.onsets)
                                        conds(mc).cond_name  = multicond.names{mc};
                                        conds(mc).onsets     = multicond.onsets{mc};
                                        conds(mc).durations  = multicond.durations{mc};
                                        conds(mc).rt_trial   = multicond.rt_trial{mc};
                                        if isfield(conds(mc),'rt_trial') && ...
                                                ~isempty(conds(mc).rt_trial)
                                            lons = length(conds(mc).onsets);
                                            lreg = length(conds(mc).rt_trial);
                                            if  lreg ~= lons
                                                out.files{1} = [];
                                                beep
                                                sprintf('Number of regression targets must be the number of trials!')
                                                disp('Please correct')
                                                return
                                            end
                                        end
                                    end
                                    design.conds = conds;
                                else
                                    design.conds = job.group(g).select.subject{j}(k).design.new_design.conds;
                                    covar = [];
                                end                           
                                ncond = length(design.conds);
                                for c = 1:ncond
                                    lons = length(design.conds(c).onsets);
                                    ldur = length(design.conds(c).durations);
                                    if ldur == 1
                                        design.conds(c).durations = repmat(design.conds(c).durations, 1, lons);
                                        ldur = length(design.conds(c).durations);
                                    end
                                    if ldur ~= lons
                                        out.files{1} = [];
                                        beep
                                        sprintf('The onsets and durations of condition %d do not have the same size!', c)
                                        disp('Please correct')
                                        return
                                    end
                                    if isfield(design.conds(c),'rt_trial') && ~isempty(design.conds(c).rt_trial)
                                        lreg = length(design.conds(c).rt_trial);
                                        if lreg ~= lons
                                            out.files{1} = [];
                                            beep
                                            sprintf('Number of regression targets must be the number of trials!')
                                            disp('Please correct')
                                            return
                                        end
                                    elseif ~isfield(design.conds(c),'rt_trial')
                                        design.conds(c).rt_trial=[];
                                    end
                                end
                                checked_conds = prt_check_design(design.conds,TR,unit,job.fmri_des.hrfover,job.fmri_des.hrfdel);
                                design.conds  = checked_conds.conds;
                                design.stats  = checked_conds.stats;
                                design.TR     = checked_conds.TR;
                                design.unit   = unit;
                                design.covar  = covar;
                                maxcond       = max([design.conds(:).scans]);
                                if nscans>1 && nscans < maxcond
                                    sprintf('Design of subject %d, group %d, modality %d, exceeds time series!',j,g,k)
                                    disp('Corresponding events were discarded')                                  
                                    for l = 1:length(design.conds)
                                        ovser = find(design.conds(l).scans > nscans);
                                        inser = find(design.conds(l).scans <= nscans);
                                        design.conds(l).discardedscans = [design.conds(l).discardedscans, design.conds(l).scans(ovser)];
                                        design.conds(l).scans = design.conds(l).scans(inser);
                                        design.conds(l).blocks = design.conds(l).blocks(inser);   
                                    end
                                end
                            end
                        end
                        % Create PRT.mat modalities
                        PRT.group(g).gr_name                        = job.group(g).gr_name;
                        PRT.group(g).subject(j).subj_name           = subj_name;
                        PRT.group(g).subject(j).modality(k)         = job.group(g).select.subject{j}(k);
                        PRT.group(g).subject(j).modality(k).TR      = job.group(g).select.subject{j}(k).TR;
                        PRT.group(g).subject(j).modality(k).design  = design;
                        PRT.group(g).subject(j).modality(k).scans   = char(job.group(g).select.subject{j}(k).scans);
                   
                    end
                end
            end
            if nmod ~= length(unique(mod_names_subj));
                out.files{1} = [];
                beep;
                sprintf('Names of modalities in subject %d group %d repeated! Please correct!',j,g)
                return
            end
        end
        PRT.group(g).hrfoverlap = job.fmri_des.hrfover;
        PRT.group(g).hrfdelay   = job.fmri_des.hrfdel;
    end
end

% Save masks at the end
% -------------------------------------------------------------------------
PRT.masks  = masks;

% Save PRT.mat file
% -------------------------------------------------------------------------
disp('Saving PRT.mat.......>>')
if spm_check_version('MATLAB','7') < 0
    save(fname,'-V6','PRT');
else
    save(fname,'PRT');
end

% Review
% -------------------------------------------------------------------------
if job.review
    prt_data_review('UserData',{PRT,job.dir_name{1}});
end
    
% Function output
% -------------------------------------------------------------------------
out.files{1} = fname;
% get the group_names
for g = 1:ngroup
    out.(sprintf('gr_name%d',g)) = PRT.group(g).gr_name;
end
% get the mod_names -> use the ones from the masks!
for m = 1:numel(mod_names_uniq)
    out.(sprintf('mod_name%d',m)) = mod_names_uniq{m};
end

disp('Done')

return

% Old code to deal with covariates per trial
% if ~isempty(job.group(g).select.subject{j}(k).design.new_design.covar{1})
%     try
%         load(char(job.group(g).select.subject{j}(k).design.new_design.covar{1}));
%         if exist('R','var')
%             if size(R,1) == nscans
%                 covar = R;
%             else
%                 out.files{1} = [];
%                 beep
%                 sprintf('Number of covariates must be the number of scans! ')
%                 disp('Please correct!')
%                 return
%             end
%         else
%             out.files{1} = [];
%             beep
%             sprintf('Covariates file must contain ''R'' variable! ')
%             disp('Please correct!')
%             return
%         end
%     catch
%         beep
%         sprintf('Could not load %s file!',char(job.group(g).select.subject{j}(k).design.new_design.covar{1}))
%         out.files{1} = [];
%         return
%     end
% else
%     
% end
