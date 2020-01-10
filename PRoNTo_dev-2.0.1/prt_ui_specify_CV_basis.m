function varargout = prt_ui_specify_CV_basis(varargin)
% PRT_UI_SPECIFY_CV_BASIS M-file for prt_ui_specify_CV_basis.fig
% 
% PRT_UI_SPECIFY_CV_BASIS, by itself, creates a new PRT_UI_SPECIFY_CV_BASIS or raises the existing
% singleton*.
%
% H = PRT_UI_SPECIFY_CV_BASIS returns the handle to a new PRT_UI_SPECIFY_CV_BASIS or the handle to
% the existing singleton*.
%
% PRT_UI_SPECIFY_CV_BASIS('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_SPECIFY_CV_BASIS.M with the given input arguments.
%
% PRT_UI_SPECIFY_CV_BASIS('Property','Value',...) creates a new PRT_UI_SPECIFY_CV_BASIS or raises the
% existing singleton*.  Starting from the left, property value pairs are
% applied to the GUI before prt_ui_specify_CV_basis_OpeningFcn gets called.  An
% unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_specify_CV_basis_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2015 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_specify_CV_basis

% Last Modified by GUIDE v2.5 10-Jun-2013 11:50:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_specify_CV_basis_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_specify_CV_basis_OutputFcn, ...
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


% --- Executes just before prt_ui_specify_CV_basis is made visible.
function prt_ui_specify_CV_basis_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_specify_CV_basis (see VARARGIN)
set(handles.figure1,'Name','PRoNTo :: Specify CV')
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

a=varargin{1};
cvlist=get(a.pop_cv,'String');
set(handles.basepop,'String',cvlist(2:end))
set(handles.basepop,'Value',1)
set(handles.nfolds,'Enable','off')
set(handles.basepop,'Enable','off')
set(handles.br_load,'Enable','off')
handles.prt=a.dat;
handles.in=a.in; %inputs for call to prt_model

%get useful info in one structure
handles.legs=a.legs;
handles.cv.type='custom';
% Choose default command line output for prt_ui_specify_CV_basis
handles.output = hObject;
handles.flagdone=0;
% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_specify_CV_basis_OutputFcn(hObject, eventdata, handles) 
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

%This figure can be deleted now
if isfield(handles,'figure1') && handles.flagdone
    delete(handles.figure1)
end


% --- Executes on button press in base.
function base_Callback(hObject, eventdata, handles)
% hObject    handle to base (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of base
a=get(handles.base,'Value');
if a
    set(handles.basepop,'Enable','on')
    handles.spectype=2;
    set(handles.specnum,'Value',0)
    set(handles.load,'Value',0)
    val=get(handles.basepop,'Value');
    mach=get(handles.basepop,'String');
    handles.cv.k=0; %by default, Leave-One-Out options
    if val==0
        warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
        set(handles.basepop,'Value',1)
        val=1;
    end
    if any(strfind(mach{val},'Subject Out'))
        handles.cv.type = 'loso';
    elseif any(strfind(mach{val},'Subject per Group'))
        handles.cv.type = 'losgo';
    elseif any(strfind(mach{val},'Block'))
        handles.cv.type = 'lobo';
    elseif any(strfind(mach{val},'Run'))        %currently implemented for MCKR only
        handles.cv.type = 'loro';
    end
    if any(strfind(mach{val},'k-fold'))
        kt=prt_text_input('Title','Specify k, the number of folds');
        handles.cv.k=str2double(kt);
    end
else
    set(handles.basepop,'Enable','off')
    set(handles.basepop,'Value',1)
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in specnum.
function specnum_Callback(hObject, eventdata, handles)
% hObject    handle to specnum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of specnum
a=get(handles.specnum,'Value');
if a
    set(handles.nfolds,'Enable','on')
    handles.spectype=3;
    set(handles.load,'Value',0)
    set(handles.base,'Value',0)
else
    set(handles.nfolds,'Enable','off')
    set(handles.nfolds,'Value',0)
end

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in basepop.
function basepop_Callback(hObject, eventdata, handles)
% hObject    handle to basepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns basepop contents as cell array
%        contents{get(hObject,'Value')} returns selected item from basepop
val=get(handles.basepop,'Value');
mach=get(handles.basepop,'String');
handles.cv.k=0; %by default, Leave-One-Out options
if val==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.basepop,'Value',1)
    val=1;
end
if any(strfind(mach{val},'Subject Out'))
    handles.cv.type = 'loso';
elseif any(strfind(mach{val},'Subject per Group'))
    handles.cv.type = 'losgo';
elseif any(strfind(mach{val},'Block'))
    handles.cv.type = 'lobo';
elseif any(strfind(mach{val},'Run'))        %currently implemented for MCKR only
    handles.cv.type = 'loro';
end
if any(strfind(mach{val},'k-fold'))
    kt=prt_text_input('Title','Specify k, the number of folds');
    handles.cv.k=str2double(kt);
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function basepop_CreateFcn(hObject, eventdata, handles)
% hObject    handle to basepop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function nfolds_Callback(hObject, eventdata, handles)
% hObject    handle to nfolds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nfolds as text
%        str2double(get(hObject,'String')) returns contents of nfolds as a
%        double
handles.cv.k=str2double(get(handles.nfolds,'String'));
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function nfolds_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nfolds (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load
a=get(handles.load,'Value');
if a
    set(handles.br_load,'Enable','on')
    handles.spectype=1;
    set(handles.specnum,'Value',0)
    set(handles.base,'Value',0)
else
    set(handles.br_load,'Enable','off')
end

% --- Executes on button press in br_load.
function br_load_Callback(hObject, eventdata, handles)
% hObject    handle to br_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cvmatf=spm_select(1,'mat','Select .mat file corresponding to the custom cross-validation');
handles.cv.type='custom';
handles.cv.mat_file = cvmatf;
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in donebut.
function donebut_Callback(hObject, eventdata, handles)
% hObject    handle to donebut (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
in=handles.in;
in.cv.type=handles.cv.type;
if isfield(handles.cv,'k')
    in.cv.k=handles.cv.k;
else
    in.cv.k=0;
end
if isfield(handles.cv,'mat_file')
    in.cv.mat_file=handles.cv.mat_file;
else
    in.cv.mat_file=[];
end
[d, CV, ID]=prt_model(handles.prt, in);
handles.flagdone=1;
delete(handles.figure1)
prt_ui_custom_CV(CV,ID,in,handles.prt,handles.legs);
