function varargout = prt_ui_design(varargin)
% PRT_UI_DESIGN M-file for prt_ui_design.fig
% 
% PRT_UI_DESIGN, by itself, creates a new PRT_UI_DESIGN or raises the 
% existing singleton*.
%
% H = PRT_UI_DESIGN returns the handle to a new PRT_UI_DESIGN or the handle
% to the existing singleton*.
%
% PRT_UI_DESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_DESIGN.M with the given input arguments.
%
% PRT_UI_DESIGN('Property','Value',...) creates a new PRT_UI_DESIGN or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before prt_ui_design_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_design_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_design

% Last Modified by GUIDE v2.5 05-Sep-2011 18:37:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_design_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_design_OutputFcn, ...
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


% --- Executes just before prt_ui_design is made visible.
function prt_ui_design_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_design (see VARARGIN)

% Choose default command line output for prt_ui_design
handles.output = hObject;
%if window already exists, just put it as the current figure
Tag='DDwin';
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
    set(handles.figure1,'Name','PRoNTo :: Data and design')
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
    handles.color=color;
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
        else
            if ~isempty(find(strcmpi(get(aa(i),'Style'),{'text',...
                    'radiobutton','checkbox'})))
                set(aa(i),'BackgroundColor',color.bg1)
            elseif ~isempty(find(strcmpi(get(aa(i),'Style'),'pushbutton')))
                set(aa(i),'BackgroundColor',color.fr)
            end
        end
        set(aa(i),'FontUnits','pixel')
        xf=get(aa(i),'FontSize');
        set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
            'Units','normalized')
    end
    
    
    handles.saved=0;
    set(handles.save_data,'ForegroundColor',handles.color.high)
    set(handles.save_data,'FontWeight','bold')
    set(handles.text6,'ForegroundColor',handles.color.high)
    
    
    handles.cgr=1; %current group
    handles.cs=1;
    handles.cm=1;
    handles.cf=1;
    handles.dat=struct('dir',[],'group',[],'design',[],'masks',[]);
    handles.modlist={};
    handles.load_fsmod=0;
end


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_design wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_design_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in group_list.
function group_list_Callback(hObject, eventdata, handles)
% hObject    handle to group_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns group_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group_list
handles.cgr=get(handles.group_list,'Value');
handles.cs=1;
handles.cm=1;
handles.cf=1;
update_display_data(hObject,handles);

handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function group_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to group_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in subjects_list.
function subjects_list_Callback(hObject, eventdata, handles)
% hObject    handle to subjects_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns subjects_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        subjects_list
handles.cs=get(handles.subjects_list,'Value');
handles.cm=1;
handles.cf=1;
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function subjects_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to subjects_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in modality_list.
function modality_list_Callback(hObject, eventdata, handles)
% hObject    handle to modality_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns modality_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modality_list
handles.cm=get(handles.modality_list,'Value');
handles.cf=1;
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function modality_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modality_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in file_list.
function file_list_Callback(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns file_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        file_list
handles.cf=get(handles.file_list,'Value');
update_display_data(hObject,handles);
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function file_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to file_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_res_dir.
function br_res_dir_Callback(hObject, eventdata, handles)
% hObject    handle to br_res_dir (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.dat.dir=uigetdir(cd,'Directory to write results');
set(handles.edit1,'String',handles.dat.dir,'FontAngle','normal')
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);


function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double
handles.dat.dir=get(handles.edit1,'String');
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in gr_add.
function gr_add_Callback(hObject, eventdata, handles)
% hObject    handle to gr_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
gname=prt_text_input('Title','Enter group name');
if isnumeric(gname)
    return
end
if ~isfield(handles,'ds')
    handles.ds=cell(1);
else
    handles.ds=[handles.ds; cell(1)];
end
ngr=length(handles.ds);
handles.cgr=ngr;
handles.dat.group(ngr).gr_name=gname;
newlist=[get(handles.group_list,'String'); {gname}];
set(handles.group_list,'String',newlist);
val=get(handles.use_scans,'Value');
ren=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(ren,'Label','Rename','Callback',@rengroup);
set(handles.group_list,'UIContextMenu',ren)
handles=guidata(hObject);
if val==1
    set(handles.use_scans,'Value',0)
    set(handles.subj_add,'enable','on')
    set(handles.subj_remove,'enable','on')
end
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);


%Function called when right-clicking on the 'rename' menu
function rengroup(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.group_list,'Value');
renam=prt_text_input('Title','Rename group');
if isempty(renam)
    return
end
handles.dat.group(val).gr_name=renam;
list=get(handles.group_list,'String');
list{val}=renam;
set(handles.group_list,'String',list)
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in gr_remove.
function gr_remove_Callback(hObject, eventdata, handles)
% hObject    handle to gr_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.group_list,'String');
ngr=length(handles.ds);
if ngr==1
    nlist=[];
    handles.dat.group=[];
    handles.cs=0;
    handles.cm=0;
    handles.cf=0;
    handles=rmfield(handles,'ds');
    set(handles.mask_list,...
        'String',{'none'},...
        'Value',1);
    set(handles.text6,'ForegroundColor',handles.color.high)
    handles.modlist={};
elseif handles.cgr==1 && ngr>1
    nlist={list{2:end}};
    handles.dat.group=handles.dat.group(2:end);
    handles.ds={handles.ds{2:end}};
    handles.cgr=1;
    handles.cs=1;
    handles.cm=1;
    handles.cf=1;
elseif handles.cgr==ngr
    nlist={list{1:end-1}};
    handles.dat.group=handles.dat.group(1:end-1);
    handles.cgr=handles.cgr-1;
    handles.ds={handles.ds{1:end-1}};
    handles.cs=1;
    handles.cm=1;
    handles.cf=1;
else
    nlist=list([1:handles.cgr-1,handles.cgr+1:end]);
    handles.ds=handles.ds([1:handles.cgr-1,handles.cgr+1:end]);
    handles.dat.group=handles.dat.group([1:handles.cgr-1,handles.cgr+1:end]);
    handles.cgr=handles.cgr-1;
    handles.cs=1;
    handles.cm=1;
    handles.cf=1;
end
set(handles.group_list,'String',nlist);
val=get(handles.use_scans,'Value');
if val==1
    set(handles.use_scans,'Value',0)
    set(handles.subj_add,'enable','on')
    set(handles.subj_remove,'enable','on')
end
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
Clear_Filenumber_Callback(hObject, eventdata, handles);
guidata(hObject, handles);


function Update_Filenumber_Callback(hObject, eventdata, handles, filenum)

str = sprintf('%d files selected',filenum);
set(handles.text8,'String',str);
% Update handles structure
guidata(hObject, handles);


function Clear_Filenumber_Callback(hObject, eventdata, handles)

set(handles.text8,'String','');
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in subj_add.
function subj_add_Callback(hObject, eventdata, handles)
% hObject    handle to subj_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'ds') || length(handles.ds)<1
    beep
    disp(['Please add at least one group before adding subjects'])
    return
end
defnam=['S',num2str(length(handles.ds{handles.cgr})+1)];
sname=prt_text_input('Title','Enter subject name','UserData',defnam);
if isnumeric(sname)
    return
elseif ~isempty(strfind(lower(sname),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
if isfield(handles.dat.group(handles.cgr),'subject') && ...
        isfield(handles.dat.group(handles.cgr).subject,'subj_name')
    namlist={handles.dat.group(handles.cgr).subject(:).subj_name};
    if any(strcmpi(namlist,sname))
        disp(['This subject has already been defined, name attributed is ', ...
            'S', num2str(length(namlist)+1)]);
        sname=['S', num2str(length(namlist)+1)];
    end
end
handles.ds{handles.cgr}=[handles.ds{handles.cgr}, cell(1)];
handles.cs=length(handles.ds{handles.cgr});
if ~isfield(handles.dat.group(handles.cgr),'subject')
    handles.dat.group(handles.cgr).subject=[];
end
handles.dat.group(handles.cgr).subject(handles.cs).subj_name=sname;
newlist=[get(handles.subjects_list,'String'); {sname}];
set(handles.subjects_list,'String',newlist);
rens=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(rens,'Label','Rename','Callback',@rensubj);
set(handles.subjects_list,'UIContextMenu',rens)
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);

%Function called when right-clicking on the 'rename' menu
function rensubj(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.subjects_list,'Value');
usesc=get(handles.use_scans,'Value');
if usesc
    beep
    disp('Can not rename subjects with the "scans" option')
    return
end
renam=prt_text_input('Title','Rename subject');
if isempty(renam)
    return
end
if ~isempty(strfind(lower(renam),'scan'))
    beep
    disp('Scan(s) is a reserved name. Please correct')
    return
end
handles.dat.group(handles.cgr).subject(val).subj_name=renam;
list=get(handles.subjects_list,'String');
list{val}=renam;
set(handles.subjects_list,'String',list)
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in use_scans.
function use_scans_Callback(hObject, eventdata, handles)
% hObject    handle to use_scans (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of use_scans    
val=get(handles.use_scans,'Value');
if ~isfield(handles.dat,'group') || isempty(handles.ds)
    beep
    disp('Please add at least one group before adding scans')
    return
end
handles.ds{handles.cgr}=cell(1);
handles.cs=1;
handles.dat.group(handles.cgr).subject=[];
if val==1
    set(handles.subj_add,'enable','off')
    set(handles.subj_remove,'enable','off')
    set(handles.subjects_list,'String',{'Scans'});
    handles.dat.group(handles.cgr).subject(handles.cs).subj_name='Scans';
else
    set(handles.subj_add,'enable','on')
    set(handles.subj_remove,'enable','on')
    set(handles.subjects_list,'String',{});
end
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in subj_remove.
function subj_remove_Callback(hObject, eventdata, handles)
% hObject    handle to subj_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.subjects_list,'String');
cgr=handles.cgr;
nsubj=length(handles.ds{cgr});
if nsubj==1
    nlist=[];
    handles.dat.group(cgr).subject=[];
    handles.cs=0;
    handles.cm=0;
    handles.cf=0;
    handles.ds{cgr}={};
    set(handles.mask_list,...
        'String',{'none'},...
        'Value',1);
    set(handles.text6,'ForegroundColor',handles.color.high)
    handles.modlist={};
    handles.dat.masks=[];
elseif handles.cs==1 && nsubj>1
    nlist={list{2:end}};
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject(2:end);
    handles.ds{cgr}={handles.ds{cgr}{2:end}};
    handles.cs=1;
    handles.cm=1;
    handles.cf=1;
elseif handles.cs==nsubj
    nlist={list{1:end-1}};
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject(1:end-1);
    handles.ds{cgr}={handles.ds{cgr}{1:end-1}};
    handles.cs=handles.cs-1;
    handles.cm=1;
    handles.cf=1;
else
    nlist=list([1:handles.cgr-1,handles.cgr+1:end]);
    handles.dat.group(cgr).subject=handles.dat.group(cgr).subject([1:handles.cs-1,handles.cs+1:end]);
    handles.ds{cgr}=handles.ds{cgr}([1:handles.cs-1,handles.cs+1:end]);
    handles.cs=handles.cs-1;
    handles.cm=1;
    handles.cf=1;
end
set(handles.subjects_list,'String',nlist);
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
Clear_Filenumber_Callback(hObject, eventdata, handles);
guidata(hObject, handles);

% --- Executes on button press in mod_add.
function mod_add_Callback(hObject, eventdata, handles)
% hObject    handle to mod_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'ds') || length(handles.ds)<1 || ...
        length(handles.ds{handles.cgr})<1
    beep
    disp(['Please add at least one group and one subject before adding modalities'])
    return
end
if isfield(handles.dat.group(handles.cgr).subject(1),'modality')
    in4=handles.dat.group(handles.cgr).subject(1).modality;
else
    in4=[];
end
in1=handles.modlist;
in2=handles.dat.group(handles.cgr).subject(handles.cs);
in3=[];
in5=handles.dat;
mod=prt_data_modality('UserData',{in1,in2,in3,in4,in5});
if isnumeric(mod)
    return
end
if isempty(mod.name)
    beep
    disp('A name should be given to the modality')
    return
end
if ~isempty(handles.modlist)
    for i=1:length(handles.modlist)
        if ~any(strcmpi(handles.modlist,mod.name))
            handles.modlist=[handles.modlist, {mod.name}];
        end
    end
else
    handles.modlist=[handles.modlist; {mod.name}];
end
handles.ds{handles.cgr}{handles.cs}=[handles.ds{handles.cgr}{handles.cs}, cell(1)];
handles.cm=length(handles.ds{handles.cgr}{handles.cs});
handles.ds{handles.cgr}{handles.cs}{handles.cm}=size(mod.scans,1);
if ~isfield(handles.dat.group(handles.cgr).subject(handles.cs),'modality')
    handles.dat.group(handles.cgr).subject(handles.cs).modality=struct([]);
end
%Create modality within the dat structure, with the fields compatible with
%the batch
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).mod_name=mod.name;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).detrend=mod.detrend;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).covar=mod.covar;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).rt_subj=mod.rt_subj;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).design=mod.design;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).scans=mod.scans;
newlist=[get(handles.modality_list,'String'); {mod.name}];
set(handles.modality_list,'String',newlist);
set(handles.modality_list,'Value',length(newlist));
if ~isempty(handles.modlist)
    set(handles.mask_list,'String',handles.modlist);
    set(handles.mask_list,'Value',length(handles.modlist));
else
    set(handles.mask_list,'String',{'none'});
    set(handles.mask_list,'Value',1);
end
if ~isempty(mod.scans)
    set(handles.file_list,'String',cellstr(mod.scans));
end
handles.cf=1;
renm=uicontextmenu;
item1=uimenu(renm,'Label','Modify','Callback',@renmod);
% Update handles structure
guidata(hObject, handles);
set(handles.modality_list,'UIContextMenu',renm)
handles=guidata(hObject);
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)

% Update handles structure
filenum = size(handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).scans,1);
Update_Filenumber_Callback(hObject, eventdata, handles, filenum);
guidata(hObject, handles);

%Function called when right-clicking on the 'rename' menu
function renmod(hObject,eventdata)
handles=guidata(hObject);
val=get(handles.modality_list,'Value');
in1=handles.modlist;
in2=handles.dat.group(handles.cgr).subject(handles.cs);
in3=val;
in4=[];
in5=handles.dat;
mod=prt_data_modality('UserData',{in1,in2,in3,in4,in5});
if ~isstruct(mod)
    return
end
%update structure
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).mod_name=mod.name;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).detrend=mod.detrend;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).covar=mod.covar;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).rt_subj=mod.rt_subj;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).design=mod.design;
handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).scans=mod.scans;

%update list of modalities and remove the former name from the modlist and
%mask_list if it is not present in any other subject
list=get(handles.modality_list,'String');
if ~strcmpi(list{val},mod.name)
    flag=0;
    if ~isempty(handles.modlist)
        for i=1:length(handles.ds)  %for each group
            for j=1:length(handles.ds{i}) %for each subject
                if isfield(handles.dat.group(i).subject(j).modality,'mod_name') && ...
                    any(strcmpi({handles.dat.group(i).subject(j).modality(:).mod_name},list{val}))
                    flag=flag+1;
                end
            end
        end
        if ~flag
            for i=1:length(handles.modlist)
                if strcmpi(handles.modlist{i},list{get(handles.modality_list,'Value')})
                    if i==1
                        handles.modlist={handles.modlist{2:end}};
                    elseif i==length(handles.modlist)
                        handles.modlist={handles.modlist{1:end-1}};
                    else
                        handles.modlist={handles.modlist{1:i-1}; handles.modlist{i+1:end}};
                    end
                    break
                end
            end
        end
    end
    list{val}=mod.name;
    set(handles.modality_list,'String',list);
    %update mask list and structure
    if isempty(handles.modlist)
        set(handles.mask_list,'String',{'none'});
        set(handles.mask_list,'Value',1);
    else
        set(handles.mask_list,'String',handles.modlist);
        set(handles.mask_list,'Value',length(handles.modlist));
    end
    if isfield(handles.dat.masks,'mod_name')
        for i=1:size(handles.dat.masks,2)
            if strcmpi(handles.dat.masks(i).mod_name, list{val})
                indm=i;
            end
        end
        if indm==1
            if length(handles.dat.masks)==1
                handles.dat.masks=[];
            else
                handles.dat.masks=handles.dat.masks(2:end);
            end
        elseif indm==length(handles.dat.masks)
            handles.dat.masks=handles.dat.masks(1:end-1);
        else
            handles.dat.masks=[handles.dat.masks(1:indm-1), handles.dat.masks(indm+1:end)];
        end
    end
end
% Update handles structure
filenum = size(handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).scans,1);
Update_Filenumber_Callback(hObject, eventdata, handles, filenum);
guidata(hObject, handles);


% --- Executes on button press in mod_remove.
function mod_remove_Callback(hObject, eventdata, handles)
% hObject    handle to mod_remove (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
list=get(handles.modality_list,'String');
cgr=handles.cgr;
cs=handles.cs;
nmod=length(handles.ds{cgr}{cs});
flag=0;
if nmod==1
    nlist=[];
    handles.dat.group(cgr).subject(cs).modality=[];
    handles.cm=0;
    handles.cf=0;
    handles.ds{cgr}{cs}={};
elseif handles.cm==1 && nmod>1
    nlist={list{2:end}};
    handles.dat.group(cgr).subject(cs).modality=handles.dat.group(cgr).subject(cs).modality(2:end);
    handles.ds{cgr}{cs}={handles.ds{cgr}{cs}{2:end}};
    handles.cm=1;
    handles.cf=1;
elseif handles.cm==nmod
    nlist={list{1:end-1}};
    handles.dat.group(cgr).subject(cs).modality=handles.dat.group(cgr).subject(cs).modality(1:end-1);
    handles.ds{cgr}{cs}={handles.ds{cgr}{cs}{1:end-1}};
    handles.cm=handles.cm-1;
    handles.cf=1;
else
    nlist={list{1:handles.cm-1};list{handles.cm+1:end}};
    handles.dat.group(cgr).subject(cs).modality=handles.dat.group(cgr).subject(cs).modality([1:handles.cm-1,handles.cm+1:end]);
    handles.ds{cgr}{cs}={handles.ds{cgr}{cs}{1:handles.cm-1},handles.ds{cgr}{cs}{handles.cm+1:end}};
    handles.cm=handles.cm-1;
    handles.cf=1;
end 
set(handles.modality_list,'String',nlist);
if ~isempty(handles.modlist)
    for i=1:length(handles.ds)  %for each group
        for j=1:length(handles.ds{i}) %for each subject
            if isfield(handles.dat.group(i).subject(j).modality,'mod_name') && ...
                any(strcmpi({handles.dat.group(i).subject(j).modality(:).mod_name},list{get(handles.modality_list,'Value')}))
                flag=flag+1;
            end
        end
    end
    if ~flag
        for i=1:length(handles.modlist)
            if strcmpi(handles.modlist{i},list{get(handles.modality_list,'Value')})
                if i==1
                    handles.modlist={handles.modlist{2:end}};
                elseif i==length(handles.modlist)
                    handles.modlist={handles.modlist{1:end-1}};
                else
                    handles.modlist={handles.modlist{1:i-1}; handles.modlist{i+1:end}};
                end
                break
            end
        end
    end
end
%update mask list and structure
if isempty(handles.modlist)
    set(handles.mask_list,'String',{'none'});
    set(handles.mask_list,'Value',1);
else
    set(handles.mask_list,'String',handles.modlist);
    set(handles.mask_list,'Value',length(handles.modlist));
end
if isfield(handles.dat.masks,'mod_name')
    for i=1:size(handles.dat.masks,2)
        if strcmpi(handles.dat.masks(i).mod_name, list{get(handles.modality_list,'Value')})
            indm=i;
        end
    end
    if indm==1
        if length(handles.dat.masks)==1
            handles.dat.masks=[];
        else
            handles.dat.masks=handles.dat.masks(2:end);
        end
    elseif indm==length(handles.dat.masks)
        handles.dat.masks=handles.dat.masks(1:end-1);
    else
        handles.dat.masks=[handles.dat.masks(1:indm-1), handles.dat.masks(indm+1:end)];
    end
end
update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
if size(handles.dat.group(cgr).subject(cs).modality,1) == 0
    Clear_Filenumber_Callback(hObject, eventdata, handles);
else
    filenum = size(handles.dat.group(handles.cgr).subject(handles.cs).modality(handles.cm).scans,1);
    Update_Filenumber_Callback(hObject, eventdata, handles, filenum);
end
guidata(hObject, handles);

% --- Executes on button press in file_add.
function file_add_Callback(hObject, eventdata, handles)
% hObject    handle to file_add (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check that all previous fields were filled
if ~isfield(handles,'ds') || length(handles.ds)<1 || ...
        length(handles.ds{handles.cgr})<1 || ...
        length(handles.ds{handles.cgr}{handles.cs})<1
   beep
   disp('Please, select at least one group, subject and modality before adding files')
   return
end
cgr=handles.cgr;
cs=handles.cs;
cm=handles.cm;
try
    prevlist=cellstr(handles.dat.group(cgr).subject(cs).modality(cm).scans);
catch
    prevlist={};
end
fnames=spm_select([1 Inf],'image','Select files for the modality',prevlist);
handles.dat.group(cgr).subject(cs).modality(cm).scans=fnames;
handles.ds{cgr}{cs}{cm}=size(fnames,1);
handles.cf=1;
set(handles.file_list,'String',cellstr(fnames));
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
filenum = size(handles.dat.group(cgr).subject(cs).modality(cm).scans,1);
Update_Filenumber_Callback(hObject, eventdata, handles, filenum);
guidata(hObject, handles);



% --- Executes on selection change in mask_list.
function mask_list_Callback(hObject, eventdata, handles)
% hObject    handle to mask_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mask_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mask_list

warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
list=get(handles.mask_list,'String');
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if length(list)==1
    set(handles.mask_list,'Value',1)
end

val=list{get(handles.mask_list,'Value')};
flag=0;
if ~isfield(handles.dat.masks,'mod_name')
    sel='';
    indm=1;
else
    for i=1:size(handles.dat.masks,2)
        if strcmpi(handles.dat.masks(i).mod_name, val)
            sel=handles.dat.masks(i).fname;
            indm=i;
            flag=1;
        end
    end
    if ~flag
        sel='';
        indm=length(handles.dat.masks)+1;
    end
end
mname=spm_select(1,'image',['Select mask for ',val],cellstr(sel));
if ~isempty(mname)
    handles.dat.masks(indm).mod_name=val;
    handles.dat.masks(indm).fname=mname;
end

%check if the masks structure was completely filled by user
if length(handles.modlist)==length(handles.dat.masks)
    f=0;
    for i=1:length(handles.modlist)
        if ~isempty(handles.dat.masks(i))
            f=f+1;
        end
    end
    if f==length(handles.modlist)
        set(handles.text6,'ForegroundColor',[0 0 0])
    else
        set(handles.text6,'ForegroundColor',handles.color.high)
    end
end
    
handles.saved=0;
set(handles.save_data,'ForegroundColor',handles.color.high)
% Update handles structure
guidata(hObject, handles);
        

% --- Executes during object creation, after setting all properties.
function mask_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mask_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in load_butt.
function load_butt_Callback(hObject, eventdata, handles)
% hObject    handle to load_butt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%get and laod PRT.mat
prtname = spm_select(1,'mat','Select PRT.mat',[],pwd,'PRT.mat');
PRT=prt_load(prtname,1);
if isempty(PRT)
    beep
    disp('Could not load file')
    return
end
handles.dat=PRT;
set(handles.save_data,'ForegroundColor',[0 0 0])
%Get the different fields and create the handles.ds cell array as well as
%complete the fields which might be missing (previous versions) and
%initialize the handles structure

handles.saved=1;
%flag if the masks are not linked to a modality name
if ~isfield(PRT.masks,'mod_name')
    flagmask=0;
    handles.modlist={};
else
    handles.modlist={PRT.masks(:).mod_name};
    flagmask=1;
end
    
%get the groups, subjects and modalities   
if isfield(PRT,'group')
    ng=length(PRT.group);
    handles.ds=cell(ng,1);
    for i=1:ng
        if isfield(PRT.group(i),'subject')
            ns=length(PRT.group(i).subject);
            handles.ds{i}=cell(1,ns);
            for j=1:ns
                if ~isfield(PRT.group(i).subject,'subj_name')
                    handles.saved=0;
                    set(handles.save_data,'ForegroundColor',handles.color.high)
                    PRT.group(i).subject(j).subj_name=['S',num2str(j)];
                end
                if isfield(PRT.group(i).subject(j),'modality')
                    nm=length(PRT.group(i).subject(j).modality);
                    handles.ds{i}{j}=cell(1,nm);
                    for k=1:nm
                        %if the flagmask is set to 0, then complete the list of
                        %modalities
                        if ~flagmask
                            mname=PRT.group(i).subject(j).modality(k).mod_name;
                            if ~any(strcmpi(handles.modlist,mname))
                                handles.saved=0;
                                handles.modlist=[handles.modlist, {mname}];
                            end
                        end
                        handles.ds{i}{j}{k}=size(PRT.group(i).subject(j).modality(k).scans,1);
                    end
                end
            end
        end
    end
end

if ~flagmask==1
    disp('The files for masking are not linked to a modality name')
    disp('The information was erased. Please select the mask files')
    PRT.masks=struct();
    handles.saved=0;
    set(handles.save_data,'ForegroundColor',handles.color.high)
    set(handles.text6,'ForegroundColor',handles.color.high)
end
if ~isempty(handles.modlist)
    set(handles.mask_list,'String',handles.modlist);
else
    set(handles.mask_list,'String',{'none'});
end
handles.cgr=1;
handles.cs=1;
handles.cm=1;
handles.cf=1;
a=fileparts(prtname);
set(handles.edit1,'String',a)
set(handles.group_list,'String',{PRT.group(:).gr_name})

%set the 'rename' and 'modify' right-clicks
%for groups
ren=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(ren,'Label','Rename','Callback',@rengroup);
set(handles.group_list,'UIContextMenu',ren)
handles=guidata(hObject);
%for subjects
rens=uicontextmenu;
% Update handles structure
guidata(hObject, handles);
item1=uimenu(rens,'Label','Rename','Callback',@rensubj);
set(handles.subjects_list,'UIContextMenu',rens)
handles=guidata(hObject);
%for modalities
renm=uicontextmenu;
item1=uimenu(renm,'Label','Modify','Callback',@renmod);
% Update handles structure
guidata(hObject, handles);
set(handles.modality_list,'UIContextMenu',renm)
handles=guidata(hObject);

%Remove any field from previous computations
if isfield(PRT,'fs')
    beep
    disp('Fields refering to feature sets have been found')
    disp('These will be removed if modifications to the dataset are performed')
    disp('Previously computed models will also be deleted')
    disp('Be sure to change the directory if you want to keep trace of previous work')
    handles.load_fsmod=1;
else
    handles.load_fsmod=0;
end

handles.dat=PRT;
handles.dat.dir=a;
% Update handles structure
guidata(hObject, handles);

update_display_data(hObject,handles);
handles=guidata(hObject);
set(handles.text6,'ForegroundColor',[0 0 0])
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in save_data.
function save_data_Callback(hObject, eventdata, handles)
% hObject    handle to save_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.load_fsmod
    sure=prt_ui_sure;
    if ~sure || isempty(sure)
        return
    end
end

%Check that at least one group and one subject were entered
if ~isfield(handles,'ds') || length(handles.ds) <1 || ...
        length(handles.ds{1})<1 || length(handles.ds{1}{1})<1
    beep
    disp('Please, enter at least one group completed before saving')
    return
end

%Rearrange the data structure if the "Scans" option was selected for one of
%the groups
for i=1:length(handles.ds)
    if ~isempty(strfind(lower(handles.dat.group(i).subject(1).subj_name),'scan'))
        subj=struct();
        nsubj=handles.ds{i}{1}{1};
        for j=1:length(handles.dat.group(i).subject(1).modality)
            nsubj2=handles.ds{i}{1}{j};
            if nsubj ~=nsubj2
                beep
                sprintf('Number of subjects in modality %d and 1 of group %d are different ',j,i)
                disp('Please correct')
                return
            end
        end
        handles.ds{i}=cell(nsubj,1);
        for k=1:nsubj
            subj(k).subj_name=['S', num2str(k)];
            nmod=length(handles.dat.group(i).subject(1).modality);
            handles.ds{i}{k}=cell(nmod);
            for j=1:nmod
                handles.ds{i}{k}{j}=1;
                subj(k).modality(j)=handles.dat.group(i).subject(1).modality(j);
                subj(k).modality(j).scans=subj(k).modality(j).scans(k,:);
                if ~isempty(handles.dat.group(i).subject(1).modality(j).rt_subj)
                    subj(k).modality(j).rt_subj=subj(k).modality(j).rt_subj(k);
                else
                    subj(k).modality(j).rt_subj=[];
                end
                if ~isempty(handles.dat.group(i).subject(1).modality(j).covar)
                    subj(k).modality(j).covar=subj(k).modality(j).covar(k,:);
                else
                    subj(k).modality(j).covar=[];
                end
            end
        end
        handles.dat.group(i).subject=subj;
    end
end

%Check that the different groups have the same number of
%modalities and that the different subjects have the same number of
%modalities and files per modality
ng=length(handles.ds);
nm=length(handles.ds{1}{1});
ns=length(handles.ds{1});
list=get(handles.mask_list,'String');
nmask=length(list);
%check that one mask was entered for each modality
if nmask~=length(handles.dat.masks)
    beep
    sprintf('%d masks were found, while %d modalities were added',length(handles.dat.masks),nmask)
    disp('Please correct')
    return
end

for i=1:ng
    matdat=zeros(ns,nm);
    ns=length(handles.ds{i});
    for j=1:ns
        nmi=length(handles.ds{i}{1});        
        nmj=length(handles.ds{i}{j});
        if nmj~=nmi
            beep
            sprintf('Numbers of modalities in subjects 1 and %d from group %d differ', j,i)
            disp('Please correct')
            return
        elseif nmj~=nm
            beep
            sprintf('Numbers of modalities in groups 1 and %d differ \n',i)
            disp('Please correct')
            return
        elseif nmj~=nmask
            beep
            sprintf('%d modalities found for subject %d of group %d, while %d masks found \n',nmj,j,i,nmask)
            disp('Possible errors in the modalities names, please correct')
            return
        end
        for k=1:nm
            m2=find(strcmpi({handles.dat.group(i).subject(j).modality(:).mod_name},list(k)));
            matdat(j,k)=handles.ds{i}{j}{m2};
            if isstruct(handles.dat.group(i).subject(j).modality(m2).design)
                des=handles.dat.group(i).subject(j).modality(m2).design;
                maxcond=max([des.conds(:).scans]);
                if matdat(j,k)>1 && matdat(j,k)<maxcond
                    beep
                    sprintf('Design of subject %d, group %d, modality %d, exceeds time series \n',j,i,k)
                    disp('Corresponding events were discarded')
                    for l=1:length(des.conds)
                        ovser=find(des.conds(l).scans>matdat(j,k));
                        inser=find(des.conds(l).scans<=matdat(j,k));
                        des.conds(l).discardedscans=[des.conds(l).discardedscans; des.conds(l).scans(ovser)];
                        des.conds(l).scans=des.conds(l).scans(inser);
                        des.conds(l).blocks=des.conds(l).blocks(inser);
                    end
                    handles.dat.group(i).subject(j).modality(m2).design=des;
                end
            end
        end
        handles.dat.group(i).subject(j).modality(k)=handles.dat.group(i).subject(j).modality(m2);
    end
end
def=prt_get_defaults('datad');
%save the data structure
disp('Saving the data.....>>')
PRT=struct();
PRT.group=handles.dat.group;

if ~isfield(handles.dat.group(1),'hrfoverlap')
    for i=1:ng
        PRT.group(i).hrfoverlap=def.hrfw;
    end
end
if ~isfield(handles.dat.group(1),'hrfdelay')
    for i=1:ng
        PRT.group(i).hrfdelay=def.hrfd;
    end
end
PRT.masks=handles.dat.masks;

%Remove any field from previous computations if the PRT is loaded and then
%modified
if isfield(PRT,'fs')
    PRT=rmfield(PRT,'fs');
    beep
    disp('Fields refering to feature sets have been found')
    disp('These will be removed')
    disp('Be sure to change the directory if you want to keep trace of previous work')
end
if isfield(PRT,'fas')
    PRT=rmfield(PRT,'fas');
end
if isfield(PRT,'model')
    PRT=rmfield(PRT,'model');
end

if isempty(handles.dat.dir)
    beep
    disp('Please select folder to save PRT structure to')
    return
else
    resn=fullfile(handles.dat.dir,'PRT.mat');
    save(resn,'PRT')
    
    handles.saved=1;
    disp('Save Done')
    set(handles.save_data,'ForegroundColor',[0 0 0])
    
    
    cd(handles.dat.dir)
    
    % Update handles structure
    guidata(hObject, handles);
end


% --- Executes on button press in review_button.
function review_button_Callback(hObject, eventdata, handles)
% hObject    handle to review_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
try
    nm=length(handles.dat.group(1).subject(1).modality);
catch
    beep
    disp('Please enter at least one subject in one group before reviewing')
    return
end
try
    nma=length(handles.dat.masks);
catch
    beep
    disp('A mask should be specified for each modality before reviewing')
    return
end
if nm~=nma
    beep
    disp('Number of masks does not match number of modalities')
    disp('Please, correct')
    return
end
if ~(handles.saved)
    PRT=struct();
    PRT.group=handles.dat.group;
    PRT.masks=handles.dat.masks;
else
    fname=[get(handles.edit1,'String'),filesep,'PRT.mat'];
    PRT=prt_load(fname);
    if isempty(PRT)
        beep
        disp('Could not load the saved PRT.mat')
        return
    end
end
prt_data_review('UserData',{PRT,handles.dat.dir});
fname=[get(handles.edit1,'String'),filesep,'PRT.mat'];
PRT=prt_load(fname);
handles.dat.group=PRT.group;
handles.dat.masks=PRT.masks;
% Update handles structure
guidata(hObject, handles);

           

% --- Executes on button press in quit_data.
function quit_data_Callback(hObject, eventdata, handles)
% hObject    handle to quit_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if no valuable information was entered, then exit
if ~isfield(handles,'ds') || length(handles.ds) <1 || ...
        length(handles.ds{1})<1 || length(handles.ds{1}{1})<1
    % The figure can be deleted now
    delete(handles.figure1);
else
    if ~handles.saved
        beep
        disp('Modifications performed since last saving, please save again')
        return
    else
        % The figure can be deleted now
        delete(handles.figure1);
    end
end


% -------------------------------------------------------------------------
% --------------------- Subfunctions --------------------------------------
% -------------------------------------------------------------------------

function update_display_data(hObject,hand)

ng=hand.cgr;
ns=hand.cs;
nm=hand.cm;
nf=hand.cf;
%update group display
set(hand.group_list,'Value',ng)
%update subjects display
try
    set(hand.subjects_list,'String',{hand.dat.group(ng).subject(:).subj_name})
    set(hand.subjects_list,'Value',ns)
catch
    set(hand.subjects_list,'String',{})
end
%update modality display
try
    set(hand.modality_list,'String',{hand.dat.group(ng).subject(ns).modality(:).mod_name})
    set(hand.modality_list,'Value',nm)
catch
    set(hand.modality_list,'String',{})
end
%update file display
try
    set(hand.file_list,'String',{hand.dat.group(ng).subject(ns).modality(nm).scans})
    set(hand.file_list,'Value',nf)
catch
    set(hand.file_list,'String',{})
end

% Update handles structure
guidata(hObject, hand);
