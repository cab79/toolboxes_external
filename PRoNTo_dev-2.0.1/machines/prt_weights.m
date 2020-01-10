function weights = prt_weights(d,m)
% Run function to compute weights
% FORMAT weights = prt_weights(d,m)
% Inputs:
%       d   - data structure
%             (fields of .d can vary depending on weights function)
%       m   - machine structure
%           .function - function to compute weights (string)
%           .args     - function arguments
% Output:
%       weights - weights vector [Nfeatures x 1]
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M.J.Rosa and J.Mourao-Miranda
% $Id$

SANITYCHECK = true; % turn off for speed

% initial checks
%--------------------------------------------------------------------------
if SANITYCHECK == true
    if ~isempty(m)
        if isstruct(m)
            if isfield(m,'function')
                if ~ischar(m.function)
                    error('prt_weights:functionNotString',...
                        'Error: ''function'' should be a string!');
                end
                if ~exist(m.function,'file')
                    error('prt_weights:functionFileNotFound',...
                        'Error: %s function could not be found!',...
                        m.function);
                end
            else
                error('prt_weights:functionNotField',...
                    'Error: ''function'' should be a field of machine!');
            end
            if ~isfield(m,'args')
                m.args = [];
            end
        else
            error('prt_weights:machineNotStruct',...
                'Error: machine should be a structure!');
        end
    else
        error('prt_weights:machineEmpty',...
            'Error: ''machine'' cannot be empty!');
    end
    if ~isempty(d)
        if ~isstruct(d)
            error('prt_weights:dataNotStruct',...
                'Error: data should be a structure!');
        end
    else
        error('prt_weights:dataStructEmpty',...
            'Error: ''data'' struct cannot be empty!');
    end
end

% run weights
%--------------------------------------------------------------------------
fnch    = str2func(m.function);
weights = fnch(d,m.args);

% final checks
%--------------------------------------------------------------------------
if SANITYCHECK == true
    if ~iscell(weights)
        error('prt_weights:weightsNotVector',...
            'Error: weights should be a cell!');
    else
        if ~isvector(weights{1})
            error('prt_weights:weightsNotVector',...
                'Error: weights cell should contain vectors!');
        end
    end
end