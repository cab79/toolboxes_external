function varargout = prt_data_conditions(varargin)
% PRT_DATA_CONDITIONS M-file for prt_data_conditions.fig
% 
% PRT_DATA_CONDITIONS, by itself, creates a new PRT_DATA_CONDITIONS or 
% raises the existing singleton*.
%
% H = PRT_DATA_CONDITIONS returns the handle to a new PRT_DATA_CONDITIONS 
% or the handle to the existing singleton*.
%
% PRT_DATA_CONDITIONS('CALLBACK',hObject,eventData,handles,...) calls the 
% local function named CALLBACK in PRT_DATA_CONDITIONS.M with the given 
% input arguments.
%
% PRT_DATA_CONDITIONS('Property','Value',...) creates a new 
% PRT_DATA_CONDITIONS or raises the existing singleton*.  Starting from the
% left, property value pairs are applied to the GUI before 
% prt_data_conditions_OpeningFcn gets called.  An unrecognized property name
% or invalid value makes property application stop.  All inputs are passed 
% to prt_data_conditions_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Schrouff
% $Id$

% Edit the above text to modify the response to help prt_data_conditions

% Last Modified by GUIDE v2.5 30-Jan-2015 10:13:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_data_conditions_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_data_conditions_OutputFcn, ...
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


% --- Executes just before prt_data_conditions is made visible.
function prt_data_conditions_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_data_conditions (see VARARGIN)

% Choose default command line output for prt_data_conditions
handles.output = hObject;

%if window already exists, just put it as the current figure
Tag='DDcond';
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
    %build figure when it doesn't exist yet
    set(handles.figure1,'Name','PRoNTo :: Specify conditions')
    set(handles.condmenu,'String',{'Specify','From .mat file'})
    set(handles.condmenu,'Value',2)
%set size of the window, taking screen resolution and platform into account
S0= spm('WinSize','0',1);   %-Screen size (of the current monitor)
if ispc
    PF ='MS Sans Serif';
else
    PF = spm_platform('fonts');     %-Font names (for this platform)
    PF = PF.helvetica;
end
tmp  = [S0(3)/1280 (S0(4))/800];
ratio=min(tmp)*[1 1 1 1];
FS = 1 + 0.85*(min(ratio)-1);  %factor to scale the fonts
x=get(handles.figure1,'Position');
set(handles.figure1,'DefaultTextFontSize',FS*12,...
    'DefaultUicontrolFontSize',FS*12,...
    'DefaultTextFontName',PF,...
    'DefaultAxesFontName',PF,...
    'DefaultUicontrolFontName',PF)
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



handles.cond=struct();
if ~isempty(varargin) && strcmpi(varargin{1},'UserData')
    des=varargin{2}{1};
    szn=length(des.conds);    
    dat=cell(szn,4);
    for i=1:szn
        try
            dat{i,1}=des.conds(i).cond_name;
        catch
            dat{i,1}=['cond ',num2str(i)];
        end
        handles.cond(i).cond_name=dat{i,1};
        try
            temp=[];
            for j=1:length(des.conds(i).durations)
                temp=[temp, ' ',num2str(des.conds(i).durations(j),3)];
            end
            dat{i,3}=temp;
            handles.cond(i).durations=des.conds(i).durations;
        catch
            dat{i,3}='NaN';
            handles.cond(i).durations=[];
        end
        try
            temp=[];
            for j=1:length(des.conds(i).onsets)
                temp=[temp, ' ',num2str(des.conds(i).onsets(j),3)];
            end
            dat{i,2}=temp;
            handles.cond(i).onsets=des.conds(i).onsets;
        catch
            dat{i,2}='NaN';
            handles.cond(i).onsets=[];
        end
        try
            temp=[];
            for j=1:length(des.conds(i).rt_trial)
                temp=[temp, ' ',num2str(des.conds(i).rt_trial(j),3)];
            end
            dat{i,4}=temp;
            handles.cond(i).rt_trial=des.conds(i).rt_trial;
        catch
            dat{i,4}='NaN';
            handles.cond(i).rt_trial=[];
        end
    end
    set(handles.condtable,'visible','on');
    set(handles.condtable,'Data',dat);
    handles.trval=des.TR;
    set(handles.tredit,'String',num2str(des.TR));
    handles.unit=des.unit;
else
    dat={'cond1','0','0','[]'};
    set(handles.condtable,'visible','off');
    handles.trval=0;
    handles.covar=[];
    des.unit=1;
end
set(handles.condtable,'Data',dat);
set(handles.condtable,'ColumnName',{'Name','Onsets','Duration','Regression targets (trials)'});
set(handles.condtable,'ColumnEditable',[true,true,true,true]);
set(handles.condtable,'ColumnWidth',{'auto',130,130,0});
set(handles.condtable,'ColumnFormat',{'char','char','char','char'});
set(handles.pop_unit,'String',{'Seconds','Scans'});
if des.unit
    uv=1;
else
    uv=2;
end
set(handles.pop_unit,'Value',uv)
handles.unit=des.unit;

def=prt_get_defaults('datad');
if ~isempty(varargin) && length(varargin{2})>1
    PRT=varargin{2}{2};
    if isfield(PRT.group,'hrfdelay')
        handles.hrfdel=PRT.group(1).hrfdelay;
    else
        handles.hrfdel=def.hrfd;
    end
    if isfield(PRT.group,'hrfoverlap')
        handles.hrfover=PRT.group(1).hrfoverlap;
    else
        handles.hrfover=def.hrfw;
    end
else
    handles.hrfdel=def.hrfd;
    handles.hrfover=def.hrfw;
end
end
% Update handles structure
guidata(hObject, handles);

%UIWAIT makes prt_data_conditions wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_data_conditions_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if isfield(handles,'output') && ~isempty(handles.output)
    varargout{1} = handles.output;
else
    varargout{1}=[];
end

% The figure can be deleted now
if isfield(handles,'figure1')
    delete(handles.figure1);
end


% --- Executes on selection change in condmenu.
function condmenu_Callback(hObject, eventdata, handles)
% hObject    handle to condmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns condmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from condmenu
choice=get(handles.condmenu,'Value');
if choice==1
    ncond=str2double(prt_text_input('Title','Enter number of conditions'));
    if isnan(ncond)
        return
    end
    dat=cell(ncond,4);
    for i=1:ncond
        dat{i,1}=['cond ',num2str(i)];
        dat{i,2}='NaN';
        dat{i,3}='NaN';
        dat{i,4}='[]';
        handles.cond(i).cond_name=['cond ',num2str(i)];
    end
    set(handles.condtable,'visible','on');
    set(handles.condtable,'Data',dat);
else
    des = spm_select(1,'.mat','Select multiple conditions file');
    try
        load(des);
    catch
        beep
        disp('Could not load file')
        return
    end
    try
        na=names;
    catch
        beep
        disp('No "names" found in the .mat file, please select another file')
        return
    end
    try
        dur=durations;
    catch
        beep
        disp('No "durations" found in the .mat file, please select another file')
        return
    end
    try
        ons=onsets;
    catch
        beep
        disp('No "onsets" found in the .mat file, please select another file')
        return
    end
%     try
%         rt=rt_trial;
%     catch
%         disp('No regression target (rt_trial) found in the .mat file')
%         disp('Only classification techniques will be used')
%     end
    szn=length(names);    
    dat=cell(szn,4);
    for i=1:szn
        try
            dat{i,1}=na{i};
        catch
            dat{i,1}=['cond ',num2str(i)];
        end
        handles.cond(i).cond_name=dat{i,1};
        try
            dat{i,3}=num2str(dur{i},3);
            handles.cond(i).durations=dur{i};
        catch
            dat{i,3}='NaN';
            handles.cond(i).durations=[];
        end
        try
            dat{i,2}=num2str(ons{i},3);
            handles.cond(i).onsets=ons{i};
        catch
            dat{i,2}='NaN';
            handles.cond(i).onsets=[];
        end
        try
            dat{i,4}=num2str(rt_trial{i},3);
            handles.cond(i).rt_trial=rt_trial{i};
        catch
            dat{i,4}='NaN';
            handles.cond(i).rt_trial=[];
        end
    end
    set(handles.condtable,'visible','on');
    set(handles.condtable,'Data',dat);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function condmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to condmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function tredit_Callback(hObject, eventdata, handles)
% hObject    handle to tredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tredit as text
%        str2double(get(hObject,'String')) returns contents of tredit as a double
val=get(handles.tredit,'String');
eval(['handles.trval=',val,';']);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function tredit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tredit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_unit.
function pop_unit_Callback(hObject, eventdata, handles)
% hObject    handle to pop_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_unit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_unit
un=get(handles.pop_unit,'Value');
if un==0
    set(handles.pop_unit,'Value',1)
    un=1;
end
if un==1
    handles.unit=1;
else
    handles.unit=0;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_unit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_unit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in condtable.
function condtable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to condtable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)

ind=eventdata.Indices;
if ind(2)>1
    dat=eventdata.EditData;
    eval(['vect=[',dat,'];']);
    % vect is a vector - need to compute a scalar to be able to use ||
    %if isnan(vect) || vect>10^6 || ~any(size(vect)==1)
    if any(isnan(vect)) || any(vect>10^6) || ~any(size(vect)==1)
        beep
        disp('Bad formatting of values found!')
        sprintf('Please review and correct condition %d, column %d', ind(1), ind(2))
        return
    end
end
if ind(2)==2  %check the unicity of the onsets
    vectb=unique(vect);
    if length(vect)~=length(vectb)
        beep
        disp('Duplicated values found in the onsets of condition %d', ind(1))
        disp('Please correct')
        return
    end
end
if ind(2)==1
    handles.cond(ind(1)).cond_name=eventdata.EditData;
elseif ind(2)==2
    handles.cond(ind(1)).onsets=vect;
elseif ind(2)==3
    handles.cond(ind(1)).durations=vect;
elseif ind(2)==4
     % will never trigger for this version of ProNTo
    handles.cond(ind(1)).rt_trial=vect;
end
% TODO - version 2 will have support for rt_trial
if ~isfield(handles.cond(ind(1)),'rt_trial')
    handles.cond(ind(1)).rt_trial=[];
end
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
% hObject    handle to okbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check that the TR was entered
if handles.trval==0
    beep
    disp('Enter TR value before continuing')
    return
end


%get number of conditions 
ncond=length(handles.cond);

%check that for each condition, the size of the duration and onsets vectors
%are the same
for i=1:ncond
    szon=length(handles.cond(i).onsets);
    szdur=length(handles.cond(i).durations);
%     szrt=length(handles.cond(i).rt_trial);
    if szdur==1
        handles.cond(i).durations=repmat(handles.cond(i).durations, 1, szon);
        szdur=length(handles.cond(i).durations);
    end
    if szdur ~=szon
        beep
        sprintf('The onsets and durations of condition %d do not have the same size', i)
        disp('Please correct')
        return
    end
%     if szrt && szdur ~=szrt
%         beep
%         disp('The number of regression targets must be the number of trials!')
%         sprintf('Please correct condition %d',i)
%         return
%     end
end

%Check that the conditions do not overlap, either directly or when taking
%the width of the HRF into account
conds=prt_check_design(handles.cond,handles.trval,handles.unit,...
    handles.hrfover,handles.hrfdel);
conds.covar=handles.covar;
handles.output=conds;

% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1);

% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


%Deal with covariates at the events level - old code
% dat=get(handles.covedit,'String');
% eval(['vect=[',dat,'];']);
% %second option of loading: enter the name of a .mat file containing a
% %'R' variable
% if isnan(vect)
%     try
%         load(char(vect));
%     catch
%         beep
%         disp('Could not load file or read the covariate values')
%         disp('Please enter either a .mat file name or enter the values')
%         return
%     end
%     if ~exist('R','var')
%         beep
%         sprintf('Covariates file must contain ''R'' variable! ')
%         disp('Please correct!')
%         return
%     else
%         vect=R;
%     end
% end
% handles.covar=vect;

