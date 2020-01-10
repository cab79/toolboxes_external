function output = prt_machine(d,m)
% Run machine function for classification or regression
% FORMAT output = prt_machine(d,m)
% Inputs:
%   d            - structure with information about the data, with fields:
%    Mandatory fields:
%    .train      - training data (cell array of matrices of row vectors,
%                  each [Ntr x D]). each matrix contains one representation
%                  of the data. This is useful for approaches such as
%                  multiple kernel learning.
%    .test       - testing data  (cell array of matrices row vectors, each
%                  [Nte x D])
%    .tr_targets - training labels (for classification) or values (for
%                  regression) (column vector, [Ntr x 1])
%    .use_kernel - flag, is data in form of kernel matrices (true) or in 
%                  form of features (false)
%    Optional fields: the machine is respnsible for dealing with this
%                  optional fields (e.g. d.testcov)
%   m            - structure with information about the classification or
%                  regression machine to use, with fields:
%      .function - function for classification or regression (string)
%      .args     - function arguments (either a string, a matrix, or a
%                  struct). This is specific to each machine, e.g. for
%                  an L2-norm linear SVM this could be the C parameter
% Output:
%    output      - output of machine (struct).
%       Mandatory fields:
%       .predictions - predictions of classification or regression
%                      [Nte x D]
%       Optional fields: the machine is responsible for returning
%       parameters of interest. For exemple for an SVM this could be the
%       number of support vector used in the hyperplane weights computation
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M.J.Rosa and J.Richiardi
% $Id$

% TODO: make tr_targets a cell array (?)
% TODO: fix 80-cols limit in source code
% TODO: Multi-kernel learning

SANITYCHECK = true; % can turn off for "speed"

%% INPUT CHECKS
%--------------------------------------------------------------------------
fnch   = str2func(m.function);
if SANITYCHECK==true
    % Check machine struct properties
    if ~isempty(m)
        if isstruct(m)
            if isfield(m,'function')
                % TODO: This case maybe needs more cautious handling
                if ~exist(m.function,'file')
                    error('prt_machine:machineFunctionFileNotFound',...
                        ['Error: %s function could not be found!'],...
                        m.function);
                end
            else
                error('prt_machine:machineFunctionFieldNotFound',...
                    ['Error: machine structure should contain'...
                    ' ''.function'' field!']);
            end
            if ~isfield(m,'args')
                error('prt_machine:argsFieldNotFound',...
                    ['Error: machine structure should contain' ...
                    ' ''.args'' field!']);
            end
        else
            error('prt_machine:machineNotStruct',...
                'Error: machine should be a structure!');
        end
    else
        error('prt_machine:machineStructEmpty',...
            'Error: ''machine'' struct cannot be empty!');
    end
    
    %----------------------------------------------------------------------
    % Check data struct properties
    if ~isempty(d)
        % 1: BASIC: check all mandatory fields exist so we can relax later
        if ~isfield(d,'train')
            error('prt_machine:missingField_train',...
                ['Error: ''data'' struct must contain a ''train'' '...
                ' field!']);
        end
        if ~isfield(d,'test')
            error('prt_machine:missingField_test',...
                ['Error: ''data'' struct must contain a ''test'' '...
                ' field!']);
        end
        if ~isfield(d,'tr_targets')
            error('prt_machine:missingField_tr_targets',...
                ['Error: ''data'' struct must contain a ''tr_targets'' '...
                ' field!']);
        end
        if ~isfield(d,'use_kernel')
            error('prt_machine:missingField_use_kernel',...
                ['Error: ''data'' struct must contain a ''use_kernel'' '...
                ' field!']);
        end
        if ~isfield(d,'pred_type')
            error('prt_machine:missingField_pred_type',...
                ['Error: ''data'' struct must contain a ''pred_type'' '...
                ' field!']);
        end
        
        
        % 2: BASIC: check datatype of train/test sets
        if isempty(d.train) || isempty(d.test),
            error('prt_machine:TrAndTeEmpty',...
                'Error: training and testing data cannot be empty!');
        else
            if ~iscell(d.train) || ~iscell(d.test),
                error('prt_machine:TrAndTeEmpty',...
                    'Error: training and testing data should be cell arrays!');
            end
        end
        
        % 3: BASIC: check datatypes of labels
        if ~isempty(d.tr_targets)
            if isvector(d.tr_targets)
                % force targets to column vectors
                d.tr_targets   = d.tr_targets(:);
                Ntrain_lbs = length(d.tr_targets);
            else
                error('prt_machine:trainingLabelsNotVector',...
                    'Error: training labels should be a vector!');
            end
        else
            error('prt_machine:trainingLabelsEmpty',...
                'Error: training labels cannot be empty!');
        end
        
        % 4: Check data properties (over cells)
        Nk_train   = length(d.train);
        
        % 5: Check if data has more than one cell
%         if isempty(strfind(char(fnch),'MKL')) && Nk_train > 1
%             %Check that if multiple kernels, MKL was selected,
%             %otherwise add the kernels
%             tr_tmp = zeros(size(d.train{1}));
%             te_tmp = zeros(size(d.test{1}));
%             tecov_tmp = zeros(size(d.testcov{1}));
%             for j=1:Nk_train
%                 try
%                     %add kernels
%                     tp = d.train{j}; %train set
%                     tr_tmp=tr_tmp + tp;
%                     tp = d.test{j}; %test set
%                     te_tmp=te_tmp + tp;
%                     tp = d.testcov{j}; %test set covariance matrix for GP
%                     tecov_tmp=tecov_tmp + tp;
%                 catch
%                     error('prt_cv_model:KernelsWithDifferentDimensions', ...
%                         'Kernels cannot be added since they have different dimensions')
%                 end
%             end
%             d.train = {tr_tmp};
%             d.test = {te_tmp};
%             d.testcov = {tecov_tmp};
%             Nk_train = 1;
%             clear tr_tmp te_tmp tecov_tmp     
%         end
        
        %6: Check validity of machines chosen.(e.g. use SVM to do
        %regression is not valid
        if  strcmp(d.pred_type,'regression') 
            if ~any(strcmp(m.function,{'prt_machine_krr','prt_machine_rvr',...
                                       'prt_machine_gpml','prt_machine_gpr', 'prt_machine_sMKL_reg'}))
                error('prt_machine:RgressionMachineSupport',...
                    'Error: Regresion can only chose use KRR or RVR machines');
            end
        end
        
        % 7: Check datasets properties (within cells)
        for k = 1:Nk_train,
            if ~isempty(d.train{k}) && ~isempty(d.test{k})
                if (~prt_ismatrix(d.train{k}) && ~isvector(d.train{k}) ) || ...
                        (~prt_ismatrix(d.test{k}) && ~isvector(d.test{k}) )
                    error('prt_machine:TrAndTeNotMatrices',...
                        ['Error: training and testing datasets should ' ...
                        ' be either matrices or vectors!']);
                end
            else
                error('prt_machine:TrAndTeEmpty',...
                    'Error: training and testing datasest cannot be empty!');
            end
            % check dimensions
            [Ntrain Dtrain] = size(d.train{k});
            [Ntest, Dtest]  = size(d.test{k});
            % a: feature space dimension should be equal
            if ~(Dtrain==Dtest)
                error('prt_machine:DtrNotEqDte',['Error: Training and testing '...
                    'dimensions should match, but Dtrain=%d and Dtest=%d for '...
                    'dataset %d!'],Dtrain,Dtest,k);
            end
            % b: check we have as many training labels as examples
            if ~(Ntrain_lbs==Ntrain)
                error('prt_machine:NtrlbsNotEqNtr',['Error: Number of training '...
                    'examples and training labels should match, but Ntrain_lbs=%d '...
                    'and Ntrain=%d for dataset %d!'],Ntrain_lbs,Ntrain,k);
            end
            % c: if kernel check for kernel properties
            if d.use_kernel
                if ~(Ntrain==Dtrain)
                    error('prt_machine:NtrainNotEqDtrain',['Error: Training '...
                        'dimensions should match, but Ntr=%d and Dtr=%d for '...
                        'dataset %d!'],Ntrain,Dtrain,k);
                end
                if ~(Dtest==Ntrain)
                    error('prt_machine:DtestNotEqNtrain',['Error: Testing '...
                        'dimensions should match, but Dte=%d and Ntr=%d for '...
                        'dataset %d!'],Dtest,Ntrain,k);
                end    
            end
        end
    else
        error('prt_machine:dataStructEmpty',...
            'Error: data struct cannot be empty!');
    end
end % SANITYCHECK

%% Run model
%--------------------------------------------------------------------------
try
    output = fnch(d,m.args);
catch
    err = lasterror;
    err_ID=lower(err.identifier);
    err_libProblem = strfind(err_ID,'libnotfound');
    err_argsProblem = strfind(err_ID,'argsproblem');
    disp('prt_machine: machine did not run sucessfully.');
    if ~isempty(err_libProblem)
        error('prt_machine:libNotFound',['Error: the library for '...
            'machine %s could not be found on your path. '],m.function);
    elseif ~isempty(err_argsProblem)
        disp(['Error: the arguments supplied '...
            ' are invalid. ' ...
            'SOLUTION: Please follow the advice given by the machine.']);
        error('prt_machine:argsProblem',...
            'Error running machine %s: %s %s', ...
            m.function,err.identifier,err.message);
    else
        % we don't know what more to do here, pass it up
        disp(['SOLUTION: Please read the message below and attempt to' ...
            ' correct the problem, or ask the developpers for ' ...
            'assistance by copy-pasting all messages and explaining the'...
            ' exact steps that led to the problem.']);
        disp(['These kinds of issues are typically caused by Matlab '...
            'path problems.']);
        for en=numel(err.stack):-1:1
            e=err.stack(en);
            fprintf('%d : function [%s] in file [%s] at line [%d]\n',...
                en,e.name,e.file,e.line);
        end
        error('prt_machine:otherProblem',...
            'Error running machine %s: %s %s', ...
            m.function,err.identifier,err.message);
    end
end

%% OUTPUT CHECKS
%--------------------------------------------------------------------------
if SANITYCHECK==true
    
    % Check output properties
    if ~isfield(output,'predictions');
        error('prt_machine:outputNoPredictions',['Output of machine should '...
            'contain the field ''.predictions''.']);
    else
        % FIXME: multiple kernels / feature representations is unsupported
        % here
        % [afm] removed to test glm approach
        %if (size(output.predictions,1)~= Ntest)
        %    error('prt_machine:outputNpredictionsNotEqNte',['Error: Number '...
        %        'of predictions output and number of test examples should '...
        %        'match, but Npre=%d and Nte=%d !'],...
        %        size(output.predictions,1),Ntest);
        %end
    end
    
end % SANITYCHECK on output

end

%% local functions
function out = prt_ismatrix(A)
% ismatrix was not a built-in in Matlab 7.1, so do a homebrew
% implementation (based on Dan Vimont's Matlab libs at
% http://www.aos.wisc.edu/~dvimont/matlab but with short-circuit AND for
% "speed")
out=(ndims(A)==2) && (min(size(A)) ~= 1); % enable stricter check - a struct array should NOT pass.
end
