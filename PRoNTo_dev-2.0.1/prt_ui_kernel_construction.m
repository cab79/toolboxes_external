function varargout = prt_ui_kernel_construction(varargin)
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
% PRT_UI_KERNEL('Property','Value',...) creates a new PRT_UI_KERNEL or raises
% the existing singleton*.  Starting from the left, property value pairs are
% applied to the GUI before prt_ui_kernel_OpeningFcn gets called.  An
% unrecognized property name or invalid value makes property application
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

% Last Modified by GUIDE v2.5 26-Sep-2011 14:21:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_kernel_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_kernel_OutputFcn, ...
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
function prt_ui_kernel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_kernel (see VARARGIN)

% Choose default command line output for prt_ui_kernel
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_kernel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_kernel_OutputFcn(hObject, eventdata, handles) 
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
try
    load(fname)
    handles.dat=PRT;
catch
    beep
    disp('Could not load file')
    return
end
%get names of modalities
list={handles.dat.masks(:).mod_name};
set(handles.pop_mod,'String',list);
handles.modmask=cell(1,length(get(handles.pop_mod,'String')));
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
try
    load(fname)
    handles.dat=PRT;
    set(handles.edit_prt,'String',fname);
catch
    beep
    disp('Could not load file')
    return
end
%get names of modalities
list={handles.dat.masks(:).mod_name};
set(handles.pop_mod,'String',list);
handles.modmask=cell(1,length(get(handles.pop_mod,'String')));
% Update handles structure
guidata(hObject, handles);


function edit_kname_Callback(hObject, eventdata, handles)
% hObject    handle to edit_kname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_kname as text
%        str2double(get(hObject,'String')) returns contents of edit_kname as a double


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


% --- Executes on selection change in pop_build.
function pop_build_Callback(hObject, eventdata, handles)
% hObject    handle to pop_build (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_build contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_build


% --- Executes during object creation, after setting all properties.
function pop_build_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_build (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_mod.
function pop_mod_Callback(hObject, eventdata, handles)
% hObject    handle to pop_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_mod


% --- Executes during object creation, after setting all properties.
function pop_mod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_mod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_cond.
function pop_cond_Callback(hObject, eventdata, handles)
% hObject    handle to pop_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_cond contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_cond


% --- Executes during object creation, after setting all properties.
function pop_cond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in normbutt.
function normbutt_Callback(hObject, eventdata, handles)
% hObject    handle to normbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of normbutt


% --- Executes on button press in kdetrend.
function kdetrend_Callback(hObject, eventdata, handles)
% hObject    handle to kdetrend (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kdetrend

function edit_mask_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mask as text
%        str2double(get(hObject,'String')) returns contents of edit_mask as a double
fname=get(handles.mask_mod,'String');
val=get(handles.pop_mod,'Value');
handles.modmask{val}=fname;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_mask_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_mask.
function br_mask_Callback(hObject, eventdata, handles)
% hObject    handle to br_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
val=get(handles.pop_mod,'Value');
list=get(handles.pop_mod,'String');
fname=spm_select(1,'image',['Select mask for modality ',list(val)]);
handles.modmask{val}=fname;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in buildbutt.
function buildbutt_Callback(hObject, eventdata, handles)
% hObject    handle to buildbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

