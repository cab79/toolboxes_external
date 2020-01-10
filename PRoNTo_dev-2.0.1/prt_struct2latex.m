function prt_struct2latex(S)
% Function that takes in a structure S and writes down the latex code
% describing the whole structure and substructures recursively.
% The routine specifically generates the 'adv_PRTstruct.tex' file that is
% included, in the prt_manual.
%
% Bits of the code are copied/inspired by spm_latex.m from the SPM8
% distribution: http://www.fil.ion.ucl.ac.uk/spm
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by Christophe Phillips
% $Id$

if nargin<1,
    f = spm_select(1,'^PRT\.mat','Select a PRT.mat file');
    tmp = load(f);
    S = tmp.PRT;
end

fp = fopen(fullfile(prt('dir'),'manual','adv_PRTstruct.tex'),'w');
if fp==-1, sts = false; return; end;

%% Heading part
% fprintf(fp,'%% $Id$\n\n');
fprintf(fp,'\\chapter{%s}\n\\label{sec:%s}\n\\minitoc\n\n',...
    texify('PRT structure'),'PRTstruct');
fprintf(fp,'This is how the main {\\tt PRT} structure is organised.\n\n');

%% Deal with structure
fprintf(fp,'{\\tt PRT}\n');
h_skip = .2; % horizontal skip increment

i_nest = 1;
struct2tex(fp,S,h_skip,i_nest)

fclose(fp);

return

%==========================================================================
function struct2tex(fp,S,h_skip,i_nest)

if nargin<3, h_skip = .5; end

fieldn = fieldnames(S);
if i_nest<7
    % Begin list
    beg_txt = sprintf(['\\begin{list}{$\\bullet$}\n' ...
        '\t{\\setlength{\\labelsep}{.2cm}' ...
        '\\setlength{\\itemindent}{0cm}' ...
        '\\setlength{\\leftmargin}{%2.1fcm}}\n'],h_skip);
    fprintf(fp,'%s',beg_txt);
    
    % List of fields
    for ii=1:numel(fieldn)
        fprintf(fp,'\\item %s',texify(fieldn{ii}));
        if numel(S)>1
            fprintf(fp,'()\n');
        else
            fprintf(fp,'\n');
        end
        if isstruct(S(1).(fieldn{ii})) & ~isempty(S(1).(fieldn{ii}))
            struct2tex(fp,S(1).(fieldn{ii}),h_skip+.5,i_nest+1);
        end
    end
    
    % End list
    end_txt = '\end{list}';
    fprintf(fp,'%s\n',end_txt);
else
    % Too many nested lists for Latex to handle!
    
    % TODO: deal with this case!
end
return

%==========================================================================
function str = texify(str0)
st1  = strfind(str0,'/*');
en1  = strfind(str0,'*/');
st = [];
en = [];
for i=1:numel(st1),
    en1  = en1(en1>st1(i));
    if ~isempty(en1),
        st  = [st st1(i)];
        en  = [en en1(1)];
        en1 = en1(2:end);
    end;
end;

str = [];
pen = 1;
for i=1:numel(st),
    str = [str clean_latex(str0(pen:st(i)-1)) str0(st(i)+2:en(i)-1)];
    pen = en(i)+2;
end;
str = [str clean_latex(str0(pen:numel(str0)))];
return;

%==========================================================================
function str = clean_latex(str)
str  = strrep(str,'$','\$');
str  = strrep(str,'&','\&');
str  = strrep(str,'^','\^');
str  = strrep(str,'_','\_');
str  = strrep(str,'#','\#');
%str  = strrep(str,'\','$\\$');
str  = strrep(str,'|','$|$');
str  = strrep(str,'>','$>$');
str  = strrep(str,'<','$<$');
return;

%==========================================================================
