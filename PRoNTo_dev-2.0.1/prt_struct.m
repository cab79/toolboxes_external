function [PRT,flag] = prt_struct(PRT,prtdir)
%% Function to load the PRT.mat and check its integrity regarding the 
% fields that it is supposed to contain. Updates the PRT if needed.
%
% input  : PRT structure to check
%
% output : PRT structure updated
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id: prt_struct.m 803 2013-06-11 11:46:24Z cphillip $

flag = 1; % If any essential fields are missing, change flag and error
flagch = 0; % 1 if structure is updated
% Data and Design
%--------------------------------------------------------------------------
% Group
ng = {'gr_name','subject','hrfoverlap','hrfdelay'};
np = fieldnames(PRT.group(1));
cg = ismember(ng,np);
if ~all(cg) % missing fields
    flag = 0;
%     ita = find(cg==0);  Repair
%     for i = 1:length(ita)
%         for j=1:length(PRT.group)
%             PRT.group(j).(ng{ita(i)}) = [];
%             if ita(i)==1
%                 PRT.group(j).gr_name = '';
%             elseif ita(i)==3
%                 PRT.group(j).hrfoverlap = 0;
%             elseif ita(i)==4
%                 PRT.group(j).hrfdelay = 0;
%             end
%         end
%     end
end

% Subject
ng = {'subj_name','modality'};
np = fieldnames(PRT.group(1).subject(1));
cg = ismember(ng,np);
if ~all(cg) % missing fields
    flag = 0;
%     ita = find(cg==0);
%     for i = 1:length(ita)
%         for j=1:length(PRT.group)
%             for k  = 1:length(PRT.group(j).subject)
%                 PRT.group(j).subject(k).(ng{ita(i)}) = [];
%                 if ita(i)==1
%                     PRT.group(j).subject(k).subj_name = '';
%                 end
%             end
%         end
%     end
end

%Modality
ng = {'mod_name','covar','rt_subj','design','scans'};
np = fieldnames(PRT.group(1).subject(1).modality(1));
cg = ismember(ng,np);
if ~all(cg) % missing fields
    flagch = 1;
    ita = find(cg==0);
    for i = 1:length(ita)
        for j=1:length(PRT.group)
            for k  = 1:length(PRT.group(j).subject)
                for l = 1:length(PRT.group(j).subject(k).modality)
                    PRT.group(j).subject(k).modality(l).(ng{ita(i)}) = [];
                    if ita(i)==1
                        PRT.group(j).subject(k).modality(l).mod_name = '';

                    end
                end
            end
        end
    end
end

% Masks
%--------------------------------------------------------------------------
ng = {'mod_name','fname'};
np = fieldnames(PRT.masks(1));
cg = ismember(ng,np);
if ~all(cg) % missing fields
    flagch = 1;
    ita = find(cg==0);
    for i = 1:length(ita)
        for j=1:length(PRT.masks)
            PRT.masks(j).(ng{ita(i)}) = [];
        end
    end
end

% Feature set
%--------------------------------------------------------------------------
if isfield(PRT,'fs')
    % Feature set .fs
    ng = {'fs_name','k_file','id_col_names','fas','modality','id_mat','multkernel',...
        'atlas_name','multkernelROI'};
    np = fieldnames(PRT.fs(1));
    cg = ismember(ng,np);
    if ~all(cg) % missing fields
        flagch = 1;
        ita = find(cg==0);
        for i = 1:length(ita)
            for j=1:length(PRT.fs)
                PRT.fs(j).(ng{ita(i)}) = [];
                if ita(i)== 7
                    PRT.fs(j).multkernel = 0;
                elseif ita(i) == 1
                    PRT.fs(j).fs_name = '';
                elseif ita(i) == 9
                    PRT.fs(j).multkernelROI = 0;
                end
            end
        end
    end
    
    % fs.fas: indexes in the file array
    ng = {'im','ifa'};
    np = fieldnames(PRT.fs(1).fas);
    cg = ismember(ng,np);
    if ~all(cg) % missing fields
        flagch = 1;
        ita = find(cg==0);
        for i = 1:length(ita)
            for j=1:length(PRT.fs)
                for k = 1:length(PRT.fs(j).fas)
                    PRT.fs(j).fas(k).(ng{ita(i)}) = [];
                end
            end
        end
    end
    
    % fs.modality: parameters selected for each modality in GUI
    ng = {'mod_name','detrend','param_dt','mode','idfeat_fas','normalise','igood_kerns'};
    np = fieldnames(PRT.fs(1).modality);
    cg = ismember(ng,np);
    if ~all(cg) % missing fields
        flagch = 1;
        ita = find(cg==0);
        for i = 1:length(ita)
            for j=1:length(PRT.fs)
                for k = 1:length(PRT.fs(j).modality)
                    PRT.fs(j).modality(k).(ng{ita(i)}) = [];
                    if ita(i)==1
                        PRT.fs(j).modality(k).mod_name = '';
                    elseif ita(i)==2
                        PRT.fs(j).modality(k).detrend = 0;
                    elseif ita(i)==4
                        PRT.fs(j).modality(k).mode = 'all_scans';
                    elseif ita(i) == 6
                        PRT.fs(j).modality(k).normalise.type = 0;
                        PRT.fs(j).modality(k).normalise.scaling = [];                    
                    elseif ita(i)==7
                        PRT.fs(j).modality(k).igood_kerns = 1; % version 2.0
                    end
                end
            end
        end
    end

end


% File array
%--------------------------------------------------------------------------
if isfield(PRT,'fas')
    ng = {'mod_name','dat','detrend','param_dt','hdr','idfeat_img'};
    np = fieldnames(PRT.fas(1));
    cg = ismember(ng,np);
    if ~all(cg) % missing fields
        flagch = 1;
        ita = find(cg==0);
        for i = 1:length(ita)
            for j=1:length(PRT.fas)
                PRT.fas(j).(ng{ita(i)}) = [];
                if ita(i)== 3
                    PRT.fas(j).detrend = 0;
                elseif ita(i) == 1
                    PRT.fas(j).mod_name = '';
                end
            end
        end
    end   
end

% Model
%--------------------------------------------------------------------------

if isfield(PRT,'model')
    ng = {'model_name','input','output'};
    np = fieldnames(PRT.model(1));
    cg = ismember(ng,np);
    if ~all(cg) % missing fields
        flagch = 1;
        ita = find(cg==0);
        for i = 1:length(ita)
            for j=1:length(PRT.model)
                PRT.model(j).(ng{ita(i)}) = [];
                if ita(i) == 1
                    PRT.model(j).model_name = '';
                end
            end
        end
    end
    
    %model.input
    ng = {'use_kernel','type','machine','fs','samp_idx','include_allscans',...
        'targets','targ_allscans','cv_mat','cv_type','cv_k','use_nested_cv',...
        'nested_param','operations'};
    for j=1:length(PRT.model)
        if  ~isempty(PRT.model(j).input)
            np = fieldnames(PRT.model(j).input);
            cg = ismember(ng,np);
            if ~all(cg) % missing fields
                flagch = 1;
                ita = find(cg==0);
                for i = 1:length(ita)
                    for k= length(PRT.model(j).input)
                        PRT.model(j).input(k).(ng{ita(i)}) = [];
                        if ita(i) == 1
                            PRT.model(j).input(k).use_kernel = 0;
                        elseif ita(i)== 3
                            PRT.model(j).input(k).machine.function = '';
                            PRT.model(j).input(k).machine.args=[];
                        elseif ita(i) ==4
                            PRT.model(j).input(k).fs(1).fs_name = '';
                        elseif ita(i) == 6
                            PRT.model(j).input(k).include_allscans = 0;
                        elseif ita(i) == 11
                            PRT.model(j).input(k).cv_k = 0;
                        elseif ita(i) == 12
                            PRT.model(j).input(k).use_nested_cv = 0;
                        end
                    end
                end
            end
        end
    end
    
%   Dealing with model outputs
    % model.output
    ng = {'fold','stats','weight_ROI','weight_img'};
    for j=1:length(PRT.model)
        if ~isempty(PRT.model(j).output)
            np = fieldnames(PRT.model(j).output(1));
            cg = ismember(ng,np);
            if ~all(cg) % missing fields
                flagch = 1;
                ita = find(cg==0);
                for i = 1:length(ita)                    
                    for k = 1:length(PRT.model(j).output)
                        PRT.model(j).output(k).(ng{ita(i)}) = [];
                        if isfield(PRT.model(j).input,'class')
                            nclass = length(PRT.model(j).input.class);
                        else
                            nclass=1; %regression
                        end
                        winame = [prtdir,filesep,'weights_',PRT.model(j).model_name]; %potential weight image name
                        if nclass>2
                            for nc=1:nclass
                                if exist([winame,'_',num2str(nc),'.img'],'file')
                                    PRT.model(j).output(k).(ng{ita(i)})(nc) = {'weights_',PRT.model(j).model_name,'_',num2str(nc)};
                                end
                            end
                        else
                            if exist(winame,'file')
                                PRT.model(j).output(k).(ng{ita(i)}) = {'weights_',PRT.model(j).model_name};
                            end
                        end
                    end
                end
            end
        end
    end
end 

if flagch
    disp('PRT structure has been updated and saved')
end



%     
% %     ng = {'targets','predictions','stats','func_val',...
% %         'alpha','b'};
% %     np = fieldnames(PRT.model(1).output(1).fold(1));
% %     cg = ismember(ng,np);
% %     if ~all(cg) % missing fields
% %         ita = find(cg==0);
% %         for i = 1:length(ita)
% %             for j=1:length(PRT.model)
% %                 for k = 1:length(PRT.model(j).output(k))
% %                     for l = 1:length(PRT.model(j).output(k).fold(l))
% %                         PRT.model(j).output(k).fold(l).(ng{ita(i)}) = [];
% %                     end
% %                 end
% %             end
% %         end
% %     end



