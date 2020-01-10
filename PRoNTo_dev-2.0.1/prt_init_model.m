function [mid, PRT] = prt_init_model(PRT, in)
% function to initialise the model data structure
%
% FORMAT: Two modes are possible:
%     mid = prt_init_model(PRT, in)
%     [mid, PRT] = prt_init_model(PRT, in)
%
% USAGE 1:
% ------------------------------------------------------------------------
% function will return the id of a model or an error if it doesn't
% exist in PRT.mat
% Input:
% ------
% in.model_name: name of the model (string)
%
% Output:
% -------
% mid : is the identifier for the model in PRT.mat
%
% USAGE 2:
% -------------------------------------------------------------------------
% function will create the model in PRT.mat and overwrite it if it
% already exists.
%
% Input:
% ------
% in.model_name: name of the model to be created (string)
% in.use_kernel: use kernel or basis functions for this model (boolean)
% in.machine:    prediction machine to use for this model (struct)
% in.type:       'classification' or 'regression'
%
% Output:
% -------
% Populates the following fields in PRT.mat (copied from above):
% PRT.model(m).input.model_name
% PRT.model(m).input.type
% PRT.model(m).input.use_kernel
% PRT.model(m).input.machine
%
% Note: this function does not write PRT.mat. That should be done by the
%       calling function
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A. Marquand
% $Id$

% find model index
model_exists = false;
if isfield(PRT,'model')
    if any(strcmpi(in.model_name,{PRT.model(:).model_name}))
        mid = find(strcmpi(in.model_name,{PRT.model(:).model_name}));
        model_exists = true;
    else
        mid = length(PRT.model)+1;
    end
else
    mid = 1;
end

% do we want to create fields in PRT.mat?
if nargout == 1
    if model_exists
        % just display message and exit (returning id)
        disp(['Model ''',in.model_name,''' found in PRT.mat.']);
    else
        error('prt_init_model:modelNotFoundinPRT',...
            ['Model ''',in.model_name,''' not found in PRT.mat.']);
    end
else
    % initialise
    if model_exists
        warning('prt_init_model:modelAlreadyInPRT',['Model ''',in.model_name,...
            ''' already exists in PRT.mat. Overwriting...']);
    else
        disp(['Model ''',in.model_name,''' not found in PRT.mat. Creating...'])
    end
    % always overwrite the model
    PRT.model(mid).model_name       = in.model_name;
    PRT.model(mid).input.use_kernel = in.use_kernel;
    PRT.model(mid).input.type       = in.type;
    PRT.model(mid).input.machine    = in.machine;
    
    % Use nested CV to optimize hyperparameter?
    if isfield(in.cv,'nested')
        PRT.model(mid).input.use_nested_cv = in.cv.nested;
        PRT.model(mid).input.nested_param  = in.cv.nested_param;
        if ~isfield(in.cv,'type_nested')
            PRT.model(mid).input.cv_type_nested = in.cv.type;
            PRT.model(mid).input.cv_k_nested = in.cv.k;
        else
            PRT.model(mid).input.cv_type_nested = in.cv.type_nested;
            PRT.model(mid).input.cv_k_nested = in.cv.k_nested;
        end
    else
        PRT.model(mid).input.use_nested_cv = 0;
        PRT.model(mid).input.nested_param  = [];
    end
    
end

end