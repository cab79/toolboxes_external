function nt_greetings(reset)
%nt_greetings - display message the first time the toolbox is used

persistent greeted
if nargin>0; greeted=reset; return; end

if isempty(greeted)
    
    display(' ')
    display(['NoiseTools, version ',nt_version]);
    display('http://audition.ens.fr/adc/NoiseTools');
    display('Please cite relevant methods papers.');
    display(' ');
    
    greeted=1;
    
end