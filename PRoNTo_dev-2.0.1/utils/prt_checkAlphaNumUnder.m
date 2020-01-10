function out=prt_checkAlphaNumUnder(s)
% check whether a given string is alphanumerical or underscore
% FORMAT out =  prt_checkAlphaNumUnder(s)
% Inputs:
%   s - a string of arbitrary length to check
% Output:
%   out - logical 1 if the all chars in the string are alphanumerical
%         logical 0 otherwise
% 
% Based on isalpha_num in the identification toolbox
%__________________________________________________________________________
% Copyright (C) 2011 PRoNTo

%--------------------------------------------------------------------------
% Written by J.Richiardi
% $Id$

% sanity check
if ~ischar(s)
    error('Input is not a string');
end

for i=1:length(s)
    out=~isempty(strfind('ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_',upper(s(i))));
    if out==false
        break;
    end
end