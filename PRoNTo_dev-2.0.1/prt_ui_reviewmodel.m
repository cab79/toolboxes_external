function varargout = prt_ui_reviewmodel(varargin)
% PRT_UI_REVIEWMODEL M-file for prt_ui_reviewmodel.fig
% 
% PRT_UI_REVIEWMODEL, by itself, creates a new PRT_UI_REVIEWMODEL or raises
% the existing singleton*.
%
% H = PRT_UI_REVIEWMODEL returns the handle to a new PRT_UI_REVIEWMODEL or 
% the handle to the existing singleton*.
%
% PRT_UI_REVIEWMODEL('CALLBACK',hObject,eventData,handles,...) calls the 
% local function named CALLBACK in PRT_UI_REVIEWMODEL.M with the given 
% input arguments.
%
% PRT_UI_REVIEWMODEL('Property','Value',...) creates a new PRT_UI_REVIEWMODEL
% or raises the existing singleton*.  Starting from the left, property 
% value pairs are applied to the GUI before prt_ui_reviewmodel_OpeningFcn 
% gets called.  An unrecognized property name or invalid value makes 
% property application stop.  All inputs are passed to prt_ui_reviewmodel_OpeningFcn
% via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_reviewmodel

% Last Modified by GUIDE v2.5 25-Mar-2012 23:08:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_reviewmodel_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_reviewmodel_OutputFcn, ...
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


% --- Executes just before prt_ui_reviewmodel is made visible.
function prt_ui_reviewmodel_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_reviewmodel (see VARARGIN)
%if window already exists, just put it as the current figure
Tag='modelrev';
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
    
    set(handles.figure1,'Name','PRoNTo :: Review Model Specification')
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



% Choose default command line output for prt_ui_reviewmodel
handles.output = hObject;

if ~isempty(varargin{1}) && strcmpi(varargin{1},'UserData')
    handles.PRT=varargin{2}{1};
    handles.prtdir=varargin{2}{2};
else
    beep
    disp('PRT.mat needs to be provided')
    return
end

if ~isfield(handles.PRT,'model')
    beep
    disp('No model found in this PRT.mat')
    disp('No reviewing can be performed')
    delete(handles.figure1)
    return
end

listmod={handles.PRT.model(:).model_name};  
set(handles.pop_model,'String',listmod)
set(handles.pop_model,'Value',1)
%Get the indexes of the model and fs to display
indm=get(handles.pop_model,'Value');
if indm==0
    set(handles.pop_model,'Value',1)
    indm=1;
end
in.fs_name=handles.PRT.model(indm).input.fs.fs_name;
indf=prt_init_fs(handles.PRT,in);
handles.indm=indm;
handles.indf=indf;
end
% Update handles structure
guidata(hObject, handles);



% UIWAIT makes prt_ui_reviewmodel wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_reviewmodel_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
if ~isempty(handles)
    varargout{1} = handles.output;
end



% --- Executes on selection change in pop_model.
function pop_model_Callback(hObject, eventdata, handles)
% hObject    handle to pop_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_model contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_model
%Get the indexes of the model and fs to display
indm=get(handles.pop_model,'Value');
if indm==0
    set(handles.pop_model,'Value',1)
    indm=1;
end
in.fs_name=handles.PRT.model(indm).input.fs.fs_name;
indf=prt_init_fs(handles.PRT,in);
handles.indm=indm;
handles.indf=indf;
if ~isfield(handles.PRT.model(handles.indm).input,'class')
    set(handles.modelbutt,'Enable','off')
else
    set(handles.modelbutt,'Enable','on')
end
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function pop_model_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_model (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in modelbutt.
function modelbutt_Callback(hObject, eventdata, handles)
% hObject    handle to modelbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in cvbutt.
if isfield(handles.PRT.model(handles.indm).input,'class')
    prt_ui_select_class('UserData',{handles.PRT,handles.indf,handles.indm})
end


function cvbutt_Callback(hObject, eventdata, handles)
% hObject    handle to cvbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
prt_ui_reviewCV('UserData',{handles.PRT,handles.indm,...
    handles.indf,handles.prtdir})


% --- Executes on button press in kernbutt.
function kernbutt_Callback(hObject, eventdata, handles)
% hObject    handle to kernbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
kername=handles.PRT.fs(handles.indf).k_file;
try
    load([handles.prtdir,filesep,kername])
catch
    beep
    disp('Could not load kernel file')
    return
end
list = {};
for i=1:length(Phi)
    list = [list;{['Kernel ',num2str(i)]}];
end
hf=figure;
set(hf,'NumberTitle','off')
set(hf,'Name','PRoNTo :: Review kernel')
color=prt_get_defaults('color');
set(hf,'Color',color.bg1)
S0= spm('WinSize','0',1);
tmp  = [S0(3)/1280 (S0(4))/800];
ratio=min(tmp)*[1 1 1 1];
set(hf,'Resize','on')
set(hf,'Position',ratio.*[360 278 560 420])
set(hf,'Units','normalized')
c = uicontrol(hf,'Style','popupmenu',...
    'String',list,...
    'Position', ratio.*[170 365 220 50],...
    'Callback',@dispkern);
top = uicontrol(hf,'Style','text',...
    'String','Selected operations',...
    'Position', ratio.*[90 10 100 50]);
listop={'Sample averaging (within block)',...
    'Sample averaging (within subject/condition)',...
    'Mean centre features using training data',...
    'Divide data vectors by their norm',...
    'Regress out covariates (subject level)'};  %GLM for subjects only
if ~isempty(handles.PRT.model(handles.indm).input.operations)
    loper=listop(handles.PRT.model(handles.indm).input.operations);
else
    loper = {''};
end
lop = uicontrol(hf,'Style','listbox',...
    'String',loper,...
    'Position', ratio.*[190 10 200 50]);
set(c,'Units','normalized')
set(top,'Units','normalized')
set(lop,'Units','normalized')
FS = 1 + 0.85*(min(ratio)-1);  %factor to scale the fonts 
if ispc
    PF='MS Sans Serif';
else
    PF= spm_platform('fonts');     %-Font names (for this platform)
    PF=PF.helvetica;
end
set(hf,'DefaultTextFontSize',FS*12,...
        'DefaultUicontrolFontSize',FS*12,...
        'DefaultTextFontName',PF,...
        'DefaultAxesFontName',PF,...
        'DefaultUicontrolFontName',PF)
color=prt_get_defaults('color');
set(hf,'Color',color.bg1)
aa=get(hf,'children');
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uicontrol')
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
axes('Position',[0.15 0.2 0.7 0.7]);
set(c,'UserData',Phi)
imagesc(Phi{1})
colormap(jet)
colorbar

function dispkern(source,callbackdata)
ind = get(source,'Value');
Phi = get(source,'UserData');
imagesc(Phi{ind})
colormap(jet)
colorbar
