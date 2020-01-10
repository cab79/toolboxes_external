function varargout = prt_data_review(varargin)
% PRT_DATA_REVIEW M-file for prt_data_review.fig
%
% PRT_DATA_REVIEW, by itself, creates a new PRT_DATA_REVIEW or raises the 
% existing singleton*.
%
% H = PRT_DATA_REVIEW returns the handle to a new PRT_DATA_REVIEW or the 
% handle to the existing singleton*.
%
% PRT_DATA_REVIEW('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_DATA_REVIEW.M with the given input arguments.
%
% PRT_DATA_REVIEW('Property','Value',...) creates a new PRT_DATA_REVIEW or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before prt_data_review_OpeningFcn gets 
% called.  An unrecognized property name or invalid value makes property 
% application stop.  All inputs are passed to prt_data_review_OpeningFcn 
% via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Schrouff
% $Id$

% Edit the above text to modify the response to help prt_data_review

% Last Modified by GUIDE v2.5 28-Mar-2012 13:11:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_data_review_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_data_review_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before prt_data_review is made visible.
function prt_data_review_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_data_review (see VARARGIN)

% Choose default command line output for prt_data_review
handles.output = hObject;

%if window already exists, just put it as the current figure
Tag='DDrev';
F = findall(allchild(0),'Flat','Tag',Tag);
if length(F) > 1
    % Multiple Graphics windows - close all but most recent
    close(F(2:end))
    F = F(1);
    uistack(F,'top')
elseif length(F)==1
    uistack(F,'top')
else
    set(handles.figure1,'Tag',Tag)
    
set(handles.figure1,'Name','PRoNTo :: Review data and design')
%set size of the window, taking screen resolution and platform into account
S0= spm('WinSize','0',1);   %-Screen size (of the current monitor)
if ispc
    PF='MS Sans Serif';
else
    PF= spm_platform('fonts');     %-Font names (for this platform)
    PF=PF.helvetica;
end
tmp  = [S0(3)/1280 (S0(4))/800];
ratio=min(tmp)*[1 1 1 1];
FS = 1 + 0.85*(min(ratio)-1);  %factor to scale the fonts
x=get(handles.figure1,'Position');
set(handles.figure1,'Position',ratio.*x)
set(handles.figure1,'Resize','on')

% Choose the color of the different backgrounds and figure parameters
color=prt_get_defaults('color');
set(handles.figure1,'Color',color.bg1)
aa=get(handles.figure1,'children');
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        set(aa(i),'BackgroundColor',color.bg2)
        bb=get(aa(i),'children');
        if ~isempty(bb)
            for j=1:length(bb)
                if ~isempty(find(strcmpi(get(bb(j),'Style'),{'text',...
                        'radiobutton','checkbox'}))) 
                    set(bb(j),'BackgroundColor',color.bg2)
                elseif ~isempty(find(strcmpi(get(bb(j),'Style'),'pushbutton')))
                    set(bb(j),'BackgroundColor',color.fr)
                end
                set(bb(j),'FontUnits','pixel')
                xf=get(bb(j),'FontSize');
                set(bb(j),'FontSize',ceil(FS*xf),'FontName',PF,...
                    'FontUnits','normalized','Units','normalized')
            end
        end
    elseif strcmpi(get(aa(i),'type'),'uicontrol')
        if ~isempty(find(strcmpi(get(aa(i),'Style'),{'text',...
                'radiobutton','checkbox','listbox'})))
            set(aa(i),'BackgroundColor',color.bg1)
        elseif ~isempty(find(strcmpi(get(aa(i),'Style'),'pushbutton')))
            set(aa(i),'BackgroundColor',color.fr)
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        set(aa(i),'FontUnits','pixel')
        xf=get(aa(i),'FontSize');
        if ispc
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'FontUnits','normalized','Units','normalized')
        else
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'Units','normalized')
        end
    end
end


if ~isempty(varargin) && strcmpi(varargin{1},'UserData')
    %get number of groups and subjects/group
    PRT=varargin{2}{1};
    handles.prtdir=varargin{2}{2};
    ng=length(PRT.group);
    if ng==0
        beep
        disp('No group found at this point')
        disp('Please enter at least one completed group before reviewing')
        return
    end
    ns=zeros(ng,1);
    gname=cell(ng);
    handles.gname=gname;
    for i=1:ng
        if ~isfield(PRT.group(i),'subject')
            disp(['Group ' num2str(i) ' has no subject'])
            disp('Please enter at least one completed group before reviewing')
            return
        end
        ns(i)=length(PRT.group(i).subject);
        gname{i}=PRT.group(i).gr_name;
    end
    %get number of modalities and the index of those having design
    if ~isfield(PRT.group(1).subject(1),'modality')
        disp(['No modality has been defined.'])
        disp('Please enter at least one completed group before reviewing')
        return
    end
    nm=length(PRT.group(1).subject(1).modality);
    ind=[];
    list={};
    for i=1:nm
        % XXX should we check on the existence of 'design' field ?
        if isstruct(PRT.group(1).subject(1).modality(i).design)
            ind=[ind, i];
            list=[list, {PRT.group(1).subject(1).modality(i).mod_name}];
        end
    end   
else
    beep
    disp('No data structure found at this point')
    disp('Please enter at least one completed group before reviewing')
    return
end

%Set the texts of the different fields
set(handles.numgr,'String',num2str(ng))
set(handles.nummod,'String',num2str(nm));
if isempty(ind)
    set(handles.des,'String','No')
    set(handles.modlist,'String',{'None'})
    set(handles.modlist,'Visible','off')
    set(handles.condask,'Visible','off')
    set(handles.text10,'Visible','off')
    set(handles.axes2,'Visible','off')
    set(handles.axes3,'Visible','off')
    set(handles.numcond,'Visible','off')
    set(handles.pm1,'Visible','off')
    set(handles.pm2,'Visible','off')
    set(handles.intsc,'Visible','off')
    set(handles.befcor,'Visible','off')
    set(handles.aftcor,'Visible','off')
    set(handles.mbef,'Visible','off')
    set(handles.stdbef,'Visible','off')
    set(handles.maft,'Visible','off')
    set(handles.stdaft,'Visible','off')
    set(handles.hrfover_txt,'Visible','off')
    set(handles.hrfover_edit,'Visible','off')
    set(handles.txthrfdel,'Visible','off')
    set(handles.edit_hrfdel,'Visible','off')
    set(handles.uipanel4,'Visible','off')
    %Only displays the number of subjects per group when no design:
    x=get(handles.axes3,'Position');
    set(handles.axes1,'Position',[x(1),x(2)*1.5,x(3),3*x(4)]);
    x=get(handles.uipanel4,'Position');
    w=get(handles.uipanel3,'Position');
    set(handles.uipanel3,'Position',[x(1),x(2)*1.5,w(3),2.7*w(4)]);
    w=get(handles.figure1,'Position');
    set(handles.figure1,'Position',[w(1),w(2)+(2*w(4)/2.5),w(3),w(4)/2.5])
    w=get(handles.axes1,'Position');
    x=get(handles.text1,'Position');
    set(handles.text1,'Position',[x(1),(w(4)+x(4))*1.1,x(3),x(4)*3])
else
    def=prt_get_defaults('datad');
    if isfield(PRT.group,'hrfoverlap')
        set(handles.hrfover_edit,'String',num2str(PRT.group(1).hrfoverlap))
    else
        set(handles.hrfover_edit,'String',num2str(def.hrfw))
    end
    if isfield(PRT.group,'hrfdelay')
        set(handles.edit_hrfdel,'String',num2str(PRT.group(1).hrfdelay))
    else
        set(handles.edit_hrfdel,'String',num2str(def.hrfd))
    end
    set(handles.des,'String','Yes')
    set(handles.modlist,'String',list)
    set(handles.modlist,'Value',1)
    handles.ind=ind;
    set(handles.figure1,'CurrentAxes',handles.axes2)
    % Update handles structure
    guidata(hObject, handles);
    prt_disp_conditions(PRT,ind(1),handles,hObject);   
end

%Display the bar graph for the number of subjects/group
set(handles.figure1,'CurrentAxes',handles.axes1)
x=2:2:2*ng;
bar(handles.axes1,x,ns);
ylim([0 max(ns)+1])
xlim([1 max(x)+1])
set(handles.axes1,'XTickLabel',gname)
h=ylabel('Number of subjects');
set(h,'Rotation',90)
handles.PRT=PRT;
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_data_review wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_data_review_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% Get default command line output from handles structure
if isfield(handles,'output') && ~isempty(handles.output)
    varargout{1} = handles.output;
else
    varargout{1}=[];
end


% --- Closing the figure
function figure1_DeleteFcn(hObject,eventdata,handles)
% hObject    handle to datastruct (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on selection change in modlist.
function modlist_Callback(hObject, eventdata, handles)
% hObject    handle to modlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns modlist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modlist
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
list=get(handles.modlist,'String');
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if length(list)==1
    set(handles.modlist,'Value',1)
end

val=get(handles.modlist,'Value');
set(handles.figure1,'CurrentAxes',handles.axes2)
prt_disp_conditions(handles.PRT,handles.ind(val),handles,hObject);

% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function modlist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modlist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function hrfover_edit_Callback(hObject, eventdata, handles)
% hObject    handle to hrfover_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hrfover_edit as text
%        str2double(get(hObject,'String')) returns contents of hrfover_edit as a double
val=str2double(get(handles.hrfover_edit,'String'));
del=str2double(get(handles.edit_hrfdel,'String'));
list=get(handles.modlist,'String');
cm=get(handles.modlist,'Value');
PRT=handles.PRT;
for i=1:length(PRT.group)
    for j=1:length(PRT.group(i).subject)
        m=find(strcmpi({PRT.masks(:).mod_name},list{cm}));
        dess=PRT.group(i).subject(j).modality(m).design;
        desn=prt_check_design(dess.conds,dess.TR,dess.unit,val,del);
%         desn.covar=dess.covar;
        maxcond=max([desn.conds(:).scans]);
        lfiles=size(PRT.group(i).subject(j).modality(m).scans,1);
        if lfiles<maxcond
            sprintf('Design of subject %d, group %d, modality %d, exceeds time series \n',i,j,m)
            disp('Corresponding events were discarded')
            for l=1:length(desn.conds)
                ovser=find(desn.conds(l).scans>lfiles);
                inser=find(desn.conds(l).scans<lfiles);
                desn.conds(l).discardedscans=[desn.conds(l).discardedscans, desn.conds(l).scans(ovser)];
                desn.conds(l).scans=desn.conds(l).scans(inser);
                desn.conds(l).blocks=desn.conds(l).blocks(inser);
            end
        end
        PRT.group(i).subject(j).modality(m).design=desn;
    end
end
for g = 1:length(PRT.group)
    PRT.group(g).hrfoverlap=val;
    PRT.group(g).hrfdelay=del;
end
save([handles.prtdir,filesep,'PRT.mat'],'PRT')
disp('Design in PRT.mat updated')
if isfield(PRT,'fs')
    beep
    disp('Feature sets found in the PRT')
    disp('These do not correspond to the updated onsets')
    disp('Please, compute them anew')
end
handles.PRT=PRT;
% Update handles structure
guidata(hObject, handles);
set(handles.figure1,'CurrentAxes',handles.axes2)
prt_disp_conditions(PRT,handles.ind(cm),handles,hObject)
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function hrfover_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hrfover_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_hrfdel_Callback(hObject, eventdata, handles)
% hObject    handle to edit_hrfdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_hrfdel as text
%        str2double(get(hObject,'String')) returns contents of edit_hrfdel as a double
del=str2double(get(handles.edit_hrfdel,'String'));
val=str2double(get(handles.hrfover_edit,'String'));
list=get(handles.modlist,'String');
cm=get(handles.modlist,'Value');
PRT=handles.PRT;
for i=1:length(PRT.group)
    for j=1:length(PRT.group(i).subject)
        m=find(strcmpi({PRT.masks(:).mod_name},list{cm}));
        dess=PRT.group(i).subject(j).modality(m).design;
        desn=prt_check_design(dess.conds,dess.TR,dess.unit,val,del);
%         desn.covar=dess.covar;
        maxcond=max([desn.conds(:).scans]);
        lfiles=size(PRT.group(i).subject(j).modality(m).scans,1);
        if lfiles<maxcond
            sprintf('Design of subject %d, group %d, modality %d, exceeds time series \n',i,j,m)
            disp('Corresponding events were discarded')
            for l=1:length(desn.conds)
                ovser=find(desn.conds(l).scans>lfiles);
                inser=find(desn.conds(l).scans<lfiles);
                desn.conds(l).discardedscans=[desn.conds(l).discardedscans, desn.conds(l).scans(ovser)];
                desn.conds(l).scans=desn.conds(l).scans(inser);
                desn.conds(l).blocks=desn.conds(l).blocks(inser);
            end
        end
        PRT.group(i).subject(j).modality(m).design=desn;
    end
    PRT.group(i).hrfoverlap=val;
    PRT.group(i).hrfdelay=del;
end
save([handles.prtdir,filesep,'PRT.mat'],'PRT')
disp('Design in PRT.mat updated')
if isfield(PRT,'fs')
    beep
    disp('Feature sets found in the PRT')
    disp('These do not correspond to the updated onsets')
    disp('Please, compute them anew')
end
handles.PRT=PRT;
% Update handles structure
guidata(hObject, handles);
set(handles.figure1,'CurrentAxes',handles.axes2)
prt_disp_conditions(PRT,handles.ind(cm),handles,hObject)
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_hrfdel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_hrfdel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%--------------------------------------------------------------------------
%---------------------- Subfunctions --------------------------------------
%--------------------------------------------------------------------------
function prt_disp_conditions(dat,ind,handles,hObject)


mname=dat.group(1).subject(1).modality(ind).mod_name;
nco=length(dat.group(1).subject(1).modality(ind).design.conds);
set(handles.numcond,'String',num2str(nco));


meantp=zeros(length(dat.group),nco);
stdtp=zeros(length(dat.group),nco);
meantpdisc=zeros(length(dat.group),nco);
stdtpdisc=zeros(length(dat.group),nco);
mbef=zeros(length(dat.group),1);
stdbef=zeros(length(dat.group),1);
maft=zeros(length(dat.group),1);
stdaft=zeros(length(dat.group),1);
for i=1:length(dat.group)
    nsc=zeros(length(dat.group(i).subject),nco);
    ndisc=zeros(length(dat.group(i).subject),nco);
    overlbe=zeros(length(dat.group(i).subject),1);
    overlaf=zeros(length(dat.group(i).subject),1);
    for j=1:length(dat.group(i).subject)
        indm=find(strcmpi({dat.group(i).subject(j).modality(:).mod_name},mname));
        ncond=length(dat.group(i).subject(j).modality(indm).design.conds);
        des=dat.group(i).subject(j).modality(indm).design;
        for k=1:ncond
            nsc(j,k)=length(des.conds(k).scans);
            ndisc(j,k)=length(des.conds(k).discardedscans)+length(des.conds(k).hrfdiscardedscans);
        end
        overlbe(j)=des.stats.meanovl;
        overlaf(j)=des.stats.mgoodovl;
    end
    if size(nsc,1)==1
        meantp(i,:)=nsc;
        stdtp(i,:)=zeros(size(nsc));
        meantpdisc(i,:)=ndisc;
        stdtpdisc(i,:)=zeros(size(ndisc));
        mbef(i)=overlbe;
        stdbef(i)=zeros(size(overlbe));
        maft(i)=overlaf;
        stdaft(i)=zeros(size(overlaf));
    else
        meantp(i,:)=mean(nsc);
        stdtp(i,:)=std(nsc);
        meantpdisc(i,:)=mean(ndisc);
        stdtpdisc(i,:)=std(ndisc);
        mbef(i)=mean(overlbe);
        stdbef(i)=std(overlbe);
        maft(i)=mean(overlaf);
        stdaft(i)=std(overlaf);
    end
end

%plot the results into bar graphs
set(handles.figure1,'CurrentAxes',handles.axes2)
cla
ncond=size(meantp,2);
vecty=2:2+ncond-1;
y=vecty;
for i=2:size(meantp,1)
    vecty=vecty+(ncond+1);
    y=[y;vecty];
end
x=mean(y,2);
if length(x)==1
    x=2:2+ncond-1;
    xl=max(x)+1;
    meantp=[meantp; zeros(1,ncond)];
    stdtp=[stdtp; zeros(1,ncond)];
    meantpdisc=[meantpdisc; zeros(1,ncond)];
    stdtpdisc=[stdtpdisc; zeros(1,ncond)];
    vecty=2:2+ncond-1;
    y=vecty;
    for i=2:size(meantp,1)
        vecty=vecty+(ncond+1);
        y=[y;vecty];
    end
    x=mean(y,2);
else
    xl=max(x)+ncond-1;
end
scmax=max(max(meantp(:)),max(meantpdisc(:)));
bar(handles.axes2,x,meantp,0.9);
hold on
errorbar(handles.axes2,y,meantp,stdtp,'.k')
ylim([0 1.1*scmax])
xlim([1 xl])
h=ylabel('Number of selected scans');
set(h,'Rotation',90)
set(handles.axes2,'XTickLabel',handles.gname)

set(handles.figure1,'CurrentAxes',handles.axes3)
cla
bar(handles.axes3,x,meantpdisc,0.9);
hold on
errorbar(handles.axes3,y,meantpdisc,stdtpdisc,'.k')
ylim([0 1.1*scmax])
xlim([1 xl])
h=ylabel('Number of discarded scans');
set(h,'Rotation',90)
set(handles.axes3,'XTickLabel',handles.gname)
legend({dat.group(1).subject(1).modality(ind).design.conds(:).cond_name},...
    'Location','NorthEast')


%set the overlaps before and after HRF correction
set(handles.mbef,'String',num2str(mean(mbef)));
set(handles.stdbef,'String',num2str(mean(stdbef)));
set(handles.maft,'String',num2str(mean(maft)));
set(handles.stdaft,'String',num2str(mean(stdaft)));

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function savemenu_Callback(hObject, eventdata, handles)
% hObject    handle to savemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wd=cd;
cd(handles.prtdir)
[filename, pathname] = uiputfile( ...
{'*.png','Portable Network Graphics (*.png)';...
 '*.jpeg','JPEG figure (*.jpeg)';...
 '*.tiff','Compressed TIFF figure (*.tiff)';... 
 '*.fig','Matlab figure (*.fig)';...
 '*.pdf','Color PDF file (*.pdf)';...
 '*.epsc',  'Encapsulated PostScript (*.eps)'},...
 'Save figure as','data_and_design.png');
[a,b,c]=fileparts(filename);
ext=['-d',c(2:end)];

% Set the color of the different backgrounds and figure parameters to white
cf=get(handles.figure1,'Color');
set(handles.figure1,'Color',[1,1,1])
aa=get(handles.figure1,'children');
c=zeros(length(aa),3);
cb=cell(length(aa));
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        bb=get(aa(i),'children');
        cb{i}=zeros(length(bb),3);
        if ~isempty(bb)
            for j=1:length(bb)
                try
                    cb{i}(j,:)=get(bb(j),'BackgroundColor');
                    set(bb(j),'BackgroundColor',[1 1 1]);
                end
            end
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            c(i,:)=get(aa(i),'BackgroundColor');
            set(aa(i),'BackgroundColor',[1 1 1]);
        end
    end
end

print(handles.figure1,ext,[pathname,filesep,b],'-r500')

% Set the color of the different backgrounds and figure parameters to white
set(handles.figure1,'Color',cf)
aa=get(handles.figure1,'children');
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        bb=get(aa(i),'children');
        if ~isempty(bb)
            for j=1:length(bb)
                set(bb(j),'BackgroundColor',cb{i}(j,:));
            end
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            set(aa(i),'BackgroundColor',c(i,:));
        end
    end
end
cd(wd)

