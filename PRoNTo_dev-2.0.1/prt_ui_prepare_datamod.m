function varargout = prt_ui_prepare_datamod(varargin)
% PRT_UI_KERNEL_MODALITY M-file for prt_ui_kernel_modality.fig
% 
% PRT_UI_KERNEL_MODALITY, by itself, creates a new PRT_UI_KERNEL_MODALITY 
% or raises the existing singleton*.
%
% H = PRT_UI_KERNEL_MODALITY returns the handle to a new 
% PRT_UI_KERNEL_MODALITY or the handle to the existing singleton*.
%
% PRT_UI_KERNEL_MODALITY('CALLBACK',hObject,eventData,handles,...) calls 
% the local function named CALLBACK in PRT_UI_KERNEL_MODALITY.M with the 
% given input arguments.
%
% PRT_UI_KERNEL_MODALITY('Property','Value',...) creates a new 
% PRT_UI_KERNEL_MODALITY or raises the existing singleton*.  Starting from 
% the left, property value pairs are applied to the GUI before 
% prt_ui_kernel_modality_OpeningFcn gets called.  An unrecognized property 
% name or invalid value makes property application stop.  All inputs are 
% passed to prt_ui_kernel_modality_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_kernel_modality

% Last Modified by GUIDE v2.5 25-Jul-2013 11:28:40

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_prepare_datamod_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_prepare_datamod_OutputFcn, ...
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


% --- Executes just before prt_ui_kernel_modality is made visible.
function prt_ui_prepare_datamod_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_kernel_modality (see
% VARARGIN)
    
set(handles.figure1,'Name','PRoNTo :: Specify modality to include')
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


%GUI specific initialization
set(handles.par_name,'Visible','off')
set(handles.par_value,'Visible','off')
set(handles.pop_det,'String',{'No', ...
    'Polynomial','Discrete Cosine Transform'})
set(handles.pop_det,'Value',1)
set(handles.pop_norm,'String',{'No scaling', ...
    'Specify from .mat'});
set(handles.pop_norm,'Value',1)
set(handles.froi, 'Value',0)
set(handles.froi,'Enable','on')
set(handles.edit_atlas,'Visible','on')
set(handles.br_atlas,'Visible','on')
set(handles.edit_atlas,'Enable','off')
set(handles.br_atlas,'Enable','off')
if ~isempty(varargin{1}) && strcmpi(varargin{1},'UserData')
    handles.PRT=varargin{2}{1};
    handles.nmtc=varargin{2}{2};
else
    beep
    disp('Select a PRT.mat first')
    return
end
mod_n={handles.PRT.masks(:).mod_name}; %only one modality
set(handles.pop_mod,'String',mod_n)
set(handles.pop_mod,'Value',1)
% if only one modality and no design, suppress the "all conditions" option
if length(mod_n)==1 && (isempty(handles.PRT.group(1).subject(1).modality(1).design) ...
        || ~isstruct(handles.PRT.group(1).subject(1).modality(1).design))
    set(handles.pop_cond,'String',{'All scans'})
else
    set(handles.pop_cond,'String',{'All scans','All conditions'})
end
set(handles.pop_cond,'Value',1)
handles.mod=struct('mod_name',[],'mode',[],'mask',[],'detrend',[], ...
        'param_dt',[],'normalise',[],'matnorm',[]);
handles.mod.mod_name=mod_n(1);
handles.mod.mode='all_scans';
handles.mod.detrend=0;
handles.mod.normalise=0;
handles.mod.mask=[];
handles.mod.multroi=0;
handles.mod.atlasroi=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_kernel_modality wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_prepare_datamod_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
    % The figure can be deleted now
    delete(handles.figure1);
end




function edit_mask_Callback(hObject, eventdata, handles)
% hObject    handle to edit_mask (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_mask as text
%        str2double(get(hObject,'String')) returns contents of edit_mask as a double
handles.mod.mask=get(handles.edit_mask,'String');
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
handles.mod.mask=spm_select(1,'image','Select mask for the considered modality');
set(handles.edit_mask,'String',handles.mod.mask);
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in pop_cond.
function pop_cond_Callback(hObject, eventdata, handles)
% hObject    handle to pop_cond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_cond contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_cond
val=get(handles.pop_cond,'Value');
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if val==0
    set(handles.pop_cond,'Value',1)
end
val=get(handles.pop_cond,'Value');
if val==1
    handles.mod.mode='all_scans';
else
    handles.mod.mode='all_cond';
end
% Update handles structure
guidata(hObject, handles);

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


% --- Executes on selection change in pop_det.
function pop_det_Callback(hObject, eventdata, handles)
% hObject    handle to pop_det (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_det contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_det
val=get(handles.pop_det,'Value');
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if val==0
    set(handles.pop_det,'Value',1)
end
val=get(handles.pop_det,'Value')-1;
handles.mod.detrend=val;
if val==0 % No detrend
    set(handles.par_name,'Visible','off')
    set(handles.par_value,'Visible','off')
    handles.mod.param_dt=[];
elseif val==1 % Polynomial detrend
    set(handles.par_name,'Visible','on')
    set(handles.par_value,'Visible','on')
    set(handles.par_name,'String','Order')
    set(handles.par_value,'String','1')
    handles.mod.param_dt=1;
elseif val==2 % Discrete Cosine Transform
    set(handles.par_name,'Visible','on')
    set(handles.par_value,'Visible','on')
    set(handles.par_name,'String','Highpass filter cutoff (s)')
    set(handles.par_value,'String','128')
    handles.mod.param_dt=128;    
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_det_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_det (see GCBO)
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

% Hints: contents = get(hObject,'String') returns pop_mod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_mod
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
list=get(handles.pop_mod,'String');
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if length(list)==1
    set(handles.pop_mod,'Value',1)
end
val=get(handles.pop_mod,'Value');
if val==0
    set(handles.pop_mod,'Value',1)
end
val=get(handles.pop_mod,'Value');
handles.mod.mod_name=list(val);

% if only one modality and no design, suppress the "all conditions" option
im=find(strcmpi(list(val),{handles.PRT.group(1).subject(1).modality(:).mod_name}));
if isempty(handles.PRT.group(1).subject(1).modality(1).design) ...
        || ~isstruct(handles.PRT.group(1).subject(1).modality(im).design)
    set(handles.pop_cond,'String',{'All scans'})
else
    set(handles.pop_cond,'String',{'All scans','All conditions'})
end
set(handles.pop_cond,'Value',1)

% Update handles structure
guidata(hObject, handles);

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



function par_value_Callback(hObject, eventdata, handles)
% hObject    handle to par_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of par_value as text
%        str2double(get(hObject,'String')) returns contents of par_value as
%        a double
temp=get(handles.par_value,'Value');
handles.mod.param_dt=temp;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function par_value_CreateFcn(hObject, eventdata, handles)
% hObject    handle to par_value (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_norm.
function pop_norm_Callback(hObject, eventdata, handles)
% hObject    handle to pop_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_norm contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_norm
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
val=get(handles.pop_norm,'Value');
if val==0
    set(handles.pop_norm,'Value',1)
end
val=get(handles.pop_norm,'Value')-1;
handles.mod.normalise=val;
if val==0
    handles.mod.matnorm=[];
elseif val==1
    handles.mod.matnorm=spm_select(1,'mat','Select .mat file containing scaling');
    handles.mod.normalise=2;
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_norm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_norm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% Specify atlas for defining one kernel per region

% --- Executes on button press in froi.
function froi_Callback(hObject, eventdata, handles)
% hObject    handle to froi (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of froi
handles.mod.multroi=get(handles.froi,'Value');
if handles.mod.multroi
    set(handles.edit_atlas,'Enable','on')
    set(handles.br_atlas,'Enable','on')
else
    set(handles.edit_atlas,'Enable','off')
    set(handles.br_atlas,'Enable','off')
end
% Update handles structure
guidata(hObject, handles);

function edit_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to edit_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_atlas as text
%        str2double(get(hObject,'String')) returns contents of edit_atlas as a double
handles.mod.atlasroi=get(handles.edit_atlas,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_atlas_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in br_atlas.
function br_atlas_Callback(hObject, eventdata, handles)
% hObject    handle to br_atlas (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.mod.atlasroi=spm_select(1,'image','Select atlas to build one kernel per region');
set(handles.edit_atlas,'String',handles.mod.atlasroi);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in okbutt.
function okbutt_Callback(hObject, eventdata, handles)
% hObject    handle to okbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.mod.multroi ==1 && isempty(handles.mod.atlasroi)
    beep
    disp('One kernel per region selected but no atlas provided')
    disp('Please provide atlas')
    return
end

handles.output=handles.mod;
% Update handles structure
guidata(hObject, handles);

uiresume(handles.figure1)
