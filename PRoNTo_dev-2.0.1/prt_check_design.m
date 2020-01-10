function [conds] = prt_check_design(cond,tr,units,hrfoverlap,hrfdelay)
% FORMAT [conds] = prt_check_design(cond,tr,units,hrfoverlap)
%
% Check the design and discards scans which are either overlapping between
% conditions or which do not respect a minimum time interval between
% conditions (due to the width of the HRF function).
%
% INPUT
%   - cond  :   structure containing the names, durations and onsets of the
%               conditions
%   - tr    :   interscan interval (TR)
%   - units :   1 for seconds, 0 for scans
%   - hrfoverlap : value to correct for BOLD overlap (in seconds)
%   - hrfdelay   : value to correct for BOLD delay (in seconds)
%
% OUTPUT
% the same cond structure containing supplementary fields:
%   - scans :            scans retained for further classification
%   - discardedscans:    scans discarded because they overlapped between 
%                        conditions
%   - hrfdiscardedscans: scans discarded because they didn't respect the
%                        minimum time interval between conditions
%   - blocks:            represents the grouping of the stimuli (for
%                        cross-validation)
%   - stats:             struct containing the original time intervals, the
%                        time interval with only the 'good' scans, their
%                        means and standard deviation
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% get the number of conditions
ncond=length(cond);
%get defaults
def=prt_get_defaults('datad');
%suppose design is in seconds if no units specified
if nargin<3
    units=1;
end

%no correction for HRF BOLD overlap if none specified
if nargin<4
    hrfw=def.hrfw;
else
    hrfw=hrfoverlap;
end
%no correction for HRF BOLD delay if none specified
if nargin<5
    hrfd=def.hrfd;
else
    hrfd=hrfdelay;
end    

%check the level of overlapping between the different conditions
all=[];
conds=[];
for i=1:ncond
    condsc=[];
    bl=[];
    if units
    	cs=round((cond(i).onsets+hrfd)/tr+1);
        cdur=floor(cond(i).durations/tr);
    else
        cs=round(cond(i).onsets+(hrfd/tr)+1);
        cdur=round(cond(i).durations);
    end
    cdur=max(1,cdur);
    for j=1:length(cs)
        temp=cs(j):cs(j)+cdur(j)-1;
        condsc=[condsc,temp];
        bl=[bl j*ones(1,length(temp))]; 
    end
    cond(i).scans=unique(condsc);
    cond(i).blocks=bl;
    conds=[conds;i*ones(length(cond(i).scans),1)];
    all=[all; cond(i).scans'];
end
index=1:length(all);
[allc,indexc]=unique(all);
conds=conds(indexc);
interind=setdiff(index,indexc);
%Check if the same scans are implied in different conditions
if ~isempty(interind)
    intervect=all(interind);
    [allgood,indgood]=setxor(allc,intervect);
    conds=conds(indgood);
    flag=1;
else
    allgood=allc;
    flag=0;
end
if flag
    for i=1:ncond
        [dumb,indd]=setdiff(cond(i).scans,intervect);
        cond(i).blocks=cond(i).blocks(indd);
        cond(i).discardedscans=intersect(cond(i).scans,intervect);
        cond(i).scans=dumb;
    end
    disp('The same scans are implied in different conditions')
    disp('They will be discarded')
    disp('Review the data set for more information')
else
    for i=1:ncond
        cond(i).discardedscans=[];
    end
end

%Check for the overlapping between conditions: will not allow an
%overlapping smaller than hrfw seconds (the Haemodynamic Response FWHM) / by
%the TR.
thresh=ceil(hrfw/tr);
overlap=abs(diff(allgood));
changecond=diff(conds);
chan=find(changecond~=0);

alldisc=[];
for i=1:length(chan)
    intvl=overlap(chan(i));
    if ~(intvl>=thresh)
        togap=thresh-intvl;
        disc_c=allgood(chan(i)+1):allgood(chan(i)+1)+togap-1;
        interall=intersect(allgood,disc_c);
        alldisc=[alldisc;interall];
    end
end
if ~isempty(alldisc)
    for i=1:ncond
        [dumb,indd]=setdiff(cond(i).scans,alldisc);
        cond(i).blocks=cond(i).blocks(indd);
        cond(i).hrfdiscardedscans=intersect(cond(i).scans,alldisc);
        cond(i).scans=dumb;
    end
else
    for i=1:ncond
        cond(i).hrfdiscardedscans=[];
    end
end

stats=struct();
stats.overlap=overlap;
stats.goodscans=setdiff(allgood,alldisc);
stats.discscans=alldisc;
stats.meanovl=mean(overlap);
stats.stdovl=std(overlap);
ovlgood=abs(diff(stats.goodscans));
stats.mgoodovl=mean(ovlgood);
stats.sgoodovl=std(ovlgood);
stats.goodovl=ovlgood;

%output structure
conds=struct();
conds.conds=cond;
conds.stats=stats;
conds.TR=tr;
conds.unit=units;
return