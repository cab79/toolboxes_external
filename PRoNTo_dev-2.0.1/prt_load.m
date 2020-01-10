function PRT = prt_load(fname,flag)

% Function to load the PRT.mat and check its integrity regarding the 
% kernels and feature sets that it is supposed to contain. Updates the  set
% feature name if needed.
%
% input  : name of the PRT.mat, path included
%
% output : PRT structure updated
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

try
    load(fname)
catch
    beep
    disp('Could not load file')
    PRT = [];
    return
end

if nargin<2
    flag = 0;
end

% get path
prtdir = fileparts(fname);

% for each feature set, check that the corresponding .dat is present in the
% same directory and update the name of the file array if needed
if isfield(PRT,'fas') && ~isempty(PRT.fas)
    ind = [];
    for i=1:length(PRT.fas)
        % get the name of the file array
        if ~isempty(PRT.fas(i).dat)
            fa_name=PRT.fas(i).dat.fname;
            if ~ispc
                fname = strrep(fname,'\',filesep); 
            end 
            [fadir,fan,faext] = fileparts(fa_name);    
            if ~strcmpi(fadir,prtdir) % directories of PRT and feature set are different
                if ~exist(fullfile(prtdir,[fan,faext]),'file')  % no feature set found
                    beep
                    disp(['No feature set named ',fan,' found in the PRT directory'])
                    disp('Information linked to that feature set is deleted')
                    disp('Computing the weights or using non-kernel methods using that feature set won''t be permitted')
                else  % file exists but under the new directory
                    PRT.fas(i).dat.fname = fullfile(prtdir,[fan,faext]);
                    ind = [ind,i];
                end
            else
                ind = [ind,i];
            end
        elseif isempty(PRT.fas(i).dat) && ~isempty(PRT.fas(i).mod_name) %modality there in data and design, but not built in feature set
            ind = [ind,i];
        end
    end
    if isempty(ind)
        PRT = rmfield(PRT,'fas');
        PRT = rmfield(PRT,'fs');
    elseif length(ind)~=i
        % When a model comports the modality of the deleted fas, get rid of
        % the corresponding fs and model
        PRT.fas = PRT.fas(ind);
    end
end

%Check integrity of all PRT fields, for backward compatibility
if flag==1
    [PRT,flag] = prt_struct(PRT,prtdir);
    if ~flag
        error('prt_load:EssentialFieldsMissing',['Essentials fields are missing. ',...
            'This PRT cannot be used. Data and Design should be started from scratch'])
    end
end


save([prtdir,filesep,'PRT.mat'],'PRT')
