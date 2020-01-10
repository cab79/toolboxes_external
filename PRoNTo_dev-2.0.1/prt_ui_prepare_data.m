function varargout = prt_ui_prepare_data(varargin)
% PRT_UI_KERNEL MATLAB code for prt_ui_kernel.fig
% 
% PRT_UI_KERNEL, by itself, creates a new PRT_UI_KERNEL or raises the 
% existing singleton*.
%
% H = PRT_UI_KERNEL returns the handle to a new PRT_UI_KERNEL or the handle
% to the existing singleton*.
%
% PRT_UI_KERNEL('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_KERNEL.M with the given input arguments.
%
% PRT_UI_KERNEL('Property','Value',...) creates a new PRT_UI_KERNEL or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before prt_ui_kernel_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_kernel_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_kernel

% Last Modified by GUIDE v2.5 23-Jul-2013 12:53:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_prepare_data_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_prepare_data_OutputFcn, ...
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


% --- Executes just before prt_ui_kernel is made visible.
function prt_ui_prepare_data_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_kernel (see VARARGIN)

% Choose default command line output for prt_ui_kernel
handles.output = hObject;
%if window already exists, just put it as the current figure
Tag='FSwin';
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
    
set(handles.figure1,'Name','PRoNTo :: Prepare feature set')
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
    elseif strcmpi(get(aa(i),'type'),'uicontrol')
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


set(handles.sel_mod,'Enable','on')
set(handles.multkernflag,'Enable','off')
set(handles.multkernflag,'Value',0)
handles.kname=[];
handles.flag_mm = 0;
end
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_kernel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_prepare_data_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_prt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_prt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_prt as text
%        str2double(get(hObject,'String')) returns contents of edit_prt as a double
fname=get(handles.edit_prt,'String');
if exist('PRT','var')
    clear PRT
end
try
    load(fname)
    handles.dat=PRT;
catch
    beep
    disp('Could not load file')
    return
end
%if only one modality, than fill some fields automatically
n_mod=length(PRT.group(1).subject(1).modality);
handles.modnames={PRT.masks(:).mod_name};
if n_mod==1
    try
        handles.mod=prt_ui_prepare_datamod('UserData',{PRT,1});
        set(handles.num_mod,'Value',1)
        set(handles.num_mod,'String',1)
        set(handles.sel_mod,'String',{PRT.masks(1).mod_name})
    catch
        set(handles.edit_prt,'String','');
        return
    end
    set(handles.edit_kname,'ForegroundColor',handles.color.high)
else
    set(handles.text8,'ForegroundColor',handles.color.high)
end
handles.fname=fname;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_prt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_prt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_prt.
function br_prt_Callback(hObject, eventdata, handles)
% hObject    handle to br_prt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fname=spm_select(1,'.mat','Select PRT.mat',[],pwd,'PRT.mat');
if exist('PRT','var')
    clear PRT
end
try
    load(fname)
    handles.dat=PRT;
    set(handles.edit_prt,'String',fname);
catch
    beep
    disp('Could not load file')
    return
end
%if only one modality, than fill some fields automatically
n_mod=length(PRT.group(1).subject(1).modality);
handles.modnames={PRT.masks(:).mod_name};
if n_mod==1
    try
        handles.mod=prt_ui_prepare_datamod('UserData',{PRT,1});
        set(handles.num_mod,'Value',1)
        set(handles.num_mod,'String',1)
        set(handles.sel_mod,'String',{PRT.masks(1).mod_name})
    catch
        set(handles.edit_prt,'String','');
        return
    end
    set(handles.edit_kname,'ForegroundColor',handles.color.high)
else
    set(handles.text8,'ForegroundColor',handles.color.high)
end
handles.fname=fname;
% Update handles structure
guidata(hObject, handles);


function edit_kname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kname as text
%        str2double(get(hObject,'String')) returns contents of edit_kname as a double
handles.kname=deblank(get(handles.edit_kname,'String'));
if ~prt_checkAlphaNumUnder(handles.kname)
    beep
    disp('Kernel name should be entered in alphanumeric format only')
    disp('Please correct')
    set(handles.edit_kname,'ForegroundColor',[1,0,0])
    return
else
    set(handles.edit_kname,'ForegroundColor',[0 0 0])
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function edit_kname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_kname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function num_mod_Callback(hObject, eventdata, handles)
% hObject    handle to num_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_mod as text
%        str2double(get(hObject,'String')) returns contents of num_mod as a double
val=str2double(get(handles.num_mod,'String'));
set(handles.text8,'ForegroundColor',handles.color.high)
n_mod=length(handles.modnames);
if n_mod>1 && val>1
    set(handles.multkernflag,'Enable','on')
end
%handles.mod=struct();
list=[];
%initialize for all modalities
for i=1:n_mod
    handles.mod(i)=struct('mod_name',[],'mode',[],'mask',[],'detrend',[], ...
        'param_dt',[],'normalise',[],'matnorm',[],'multroi',[],'atlasroi',[]);
end
%get information for the selected modalities 
for i=1:val
    try
        tmp=prt_ui_prepare_datamod('UserData',{handles.dat,val});
    catch
        error('prt_ui_prepare_data:NoModSpecified','No modality was specified')
    end
    if isempty(tmp)
        error('prt_ui_prepare_data:EmptyModality','No modality was specified')
    end
    list=[list,tmp.mod_name];
    set(handles.sel_mod,'String',list)
    ind=find(strcmpi(handles.modnames,tmp.mod_name));
    handles.mod(ind)=tmp;
end
set(handles.text8,'ForegroundColor',handles.color.black)
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function num_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sel_mod.
function sel_mod_Callback(hObject, eventdata, handles)
% hObject    handle to sel_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sel_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_mod


% --- Executes during object creation, after setting all properties.
function sel_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in multkernflag.
function multkernflag_Callback(hObject, eventdata, handles)
% hObject    handle to multkernflag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of multkernflag
v=get(handles.multkernflag,'Value');
handles.flag_mm=v;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in buildbutt.
function buildbutt_Callback(hObject, eventdata, handles)
% hObject    handle to buildbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

input=struct('fname',[],'kname',[],'mod',[],'flag_mm',[]);
input.fname=handles.fname;
input.flag_mm = handles.flag_mm;
if isempty(handles.kname)
    beep
    disp('Enter a name for the feature set to be saved')
    set(handles.edit_kname,'ForegroundColor',handles.color.high)
    return
end
input.fs_name=handles.kname;
if ~isfield(handles,'mod') || isempty(handles.mod)
    beep
    disp('No modality was selected to build the dataset!')
    disp('Please, enter a number of modalities to use')
    set(handles.text8,'ForegroundColor',handles.color.high)
    return
end
input.mod=handles.mod;
load(input.fname);
prt_fs(PRT,input);
delete(handles.figure1)
