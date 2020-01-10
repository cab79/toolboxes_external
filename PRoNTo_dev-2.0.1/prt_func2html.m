% Script to generate the list of .m functions into html files 
% which can be browsed around with your favourite browser.
%
% Note that this script relies on the M2HTML package which is *NOT*
% distributed with PRoNTo!
% 
% For more information, please read the M2HTML tutorial and FAQ at:
%    <http://www.artefact.tk/software/matlab/m2html/>
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

cd(prt('dir'))

subdirs = {'.', ['.',filesep,'machines'], ...
                            ['.',filesep,'batch'],['.',filesep,'utils']}; 
    % list of directories to include
isubdirs = {'.svn','manual'}; % list of directories to exclude

m2html('mfiles',subdirs, 'htmldir',fullfile('manual','html_doc'), ...
        'recursive','off');
