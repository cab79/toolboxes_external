function param = prt_cv_opt_param(PRT,ID,model_id)
% Function to pass optional (advanced) parameters into the classifier. 
%
% This is primarily used for prediction methods that need to know something
% about the experimental design that is normally not accessible to ordinary
% (i.e. generic) prediction functions (e.g. task onsets or TR). Examples of
% this kind of method include multi-class classifier using kernel
% regression (MCKR) and the machine that implements nested cross-validation
%
% Inputs:
% -------
% PRT:      data structure
% ID:       id matrix for the current cross-validation fold
% model_id: which model are we working on?
%
% Outputs:
% --------
% param.id_fold:   the id matrix for this fold
% param.model_id:  id for the model being computed
% param.PRT:       PRT data structure
%
% Notes:
% --------
% The outputs (param.xxx) are provided for use by the classifier
% 
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by A Marquand 
% $Id$

param.id_fold  = ID;
param.model_id = model_id;
param.PRT      = PRT;

end