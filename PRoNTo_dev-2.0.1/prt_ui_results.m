function varargout = prt_ui_results(varargin)
% PRT_UI_RESULTS MATLAB code for prt_ui_results.fig
%
% PRT_UI_RESULTS, by itself, creates a new PRT_UI_RESULTS or raises the
% existing singleton*.
%
% H = PRT_UI_RESULTS returns the handle to a new PRT_UI_RESULTS or the
% handle to the existing singleton*.
%
% PRT_UI_RESULTS('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_RESULTS.M with the given input arguments.
%
% PRT_UI_RESULTS('Property','Value',...) creates a new PRT_UI_RESULTS or
% raises the existing singleton*.  Starting from the left, property value
% pairs are applied to the GUI before prt_ui_results_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_results_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id$

% Edit the above text to modify the response to help prt_ui_results

% Last Modified by GUIDE v2.5 13-Feb-2013 10:55:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prt_ui_results_OpeningFcn, ...
    'gui_OutputFcn',  @prt_ui_results_OutputFcn, ...
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


% --- Executes just before prt_ui_results is made visible.
function prt_ui_results_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_results (see VARARGIN)

%if window already exists, just put it as the current figure
Tag='Results';
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
    set(handles.figure1,'Name','PRoNTo :: Results')
    set(handles.figure1,'MenuBar','figure','WindowStyle','normal');
    
    %set size of the window, taking screen resolution and platform into account
    %--------------------------------------------------------------------------
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
    
    
    color=prt_get_defaults('color');
    set(handles.figure1,'Color',color.bg1)
    aa=get(handles.figure1,'children');
    for i=1:length(aa)
        if strcmpi(get(aa(i),'type'),'uipanel')
            set(aa(i),'BackgroundColor',color.bg2)
            bb=get(aa(i),'children');
            if ~isempty(bb)
                for j=1:length(bb)
                    if strcmpi(get(bb(j),'type'),'uipanel')
                        cc=get(bb(j),'children');
                        set(bb(j),'BackgroundColor',color.bg2)
                        for k=1:length(cc)
                            if strcmpi(get(cc(k),'type'),'uipanel')
                                dd=get(cc(k),'children');
                                set(cc(k),'BackgroundColor',color.bg2)
                                for l=1:length(dd)
                                    if strcmpi(get(dd(l),'type'),'uicontrol')
                                        if ~isempty(find(strcmpi(get(dd(l),'Style'),{'text',...
                                                'radiobutton','checkbox'})))
                                            set(dd(l),'BackgroundColor',color.bg2)
                                        elseif ~isempty(find(strcmpi(get(dd(l),'Style'),'pushbutton')))
                                            set(dd(l),'BackgroundColor',color.fr)
                                        end
                                    end
                                    set(dd(l),'FontUnits','pixel')
                                    xf=get(dd(l),'FontSize');
                                    if ispc
                                        set(dd(l),'FontSize',ceil(FS*xf),'FontName',PF,...
                                            'FontUnits','normalized','Units','normalized')
                                    else
                                        set(dd(l),'FontSize',ceil(FS*xf),'FontName',PF,...
                                            'Units','normalized')
                                    end
                                end
                            elseif strcmpi(get(cc(k),'type'),'uicontrol') && ...
                                    ~isempty(find(strcmpi(get(cc(k),'Style'),{'text',...
                                    'radiobutton','checkbox'})))
                                set(cc(k),'BackgroundColor',color.bg2)
                            elseif strcmpi(get(cc(k),'type'),'uicontrol')&& ...
                                    ~isempty(find(strcmpi(get(cc(k),'Style'),'pushbutton')))
                                set(cc(k),'BackgroundColor',color.fr)
                            end
                            set(cc(k),'FontUnits','pixel')
                            xf=get(cc(k),'FontSize');
                            set(cc(k),'FontSize',ceil(FS*xf),'FontName',PF,...
                                'Units','normalized')
                        end
                    elseif strcmpi(get(bb(j),'type'),'uicontrol') && ...
                            ~isempty(find(strcmpi(get(bb(j),'Style'),{'text',...
                            'radiobutton','checkbox'})))
                        set(bb(j),'BackgroundColor',color.bg2)
                    elseif strcmpi(get(bb(j),'type'),'uicontrol') && ...
                            ~isempty(find(strcmpi(get(bb(j),'Style'),'pushbutton')))
                        set(bb(j),'BackgroundColor',color.fr)
                    end
                    set(bb(j),'FontUnits','pixel')
                    xf=get(bb(j),'FontSize');
                    set(bb(j),'FontSize',ceil(FS*xf),'FontName',PF,...
                        'Units','normalized')
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
        if ~strcmpi(get(aa(i),'type'),'uimenu')
            set(aa(i),'FontUnits','pixel')
            xf=get(aa(i),'FontSize');
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'Units','normalized')
        end
    end
    
    
    % Initialize window
    % -------------------------------------------------------------------------
    
    if ~isfield(handles,'notinit')
        
        % Load PRT.mat
        PRT     = spm_select(1,'mat','Select PRT.mat',[],pwd,'PRT.mat');
        pathdir = regexprep(PRT,'PRT.mat', '');
        handles.pathdir = pathdir;
        handles.prtdir=fileparts(PRT);
        load(PRT);
        
        % Save PRT
        handles.PRT = PRT;
        
        % Flag to load new weights
        handles.noloadw = 0;
        
        % Load model names
        if ~isfield(PRT,'model')
            error('No models found in PRT.mat!')
        end
        nmodels = length(PRT.model);
        mi  = [];
        nmi = 0;
        for m = 1:nmodels
            if isfield(PRT.model(m),'input') && ~isempty(PRT.model(m).input)
                if isfield(PRT.model(m),'output') && ~isempty(PRT.model(m).output)
                    nmi = nmi +1;
                    model_name{nmi} = PRT.model(m).model_name;
                    mi = [mi, m];
                else
                    beep;
                    disp(sprintf('Model %s not estimated! It will not be displayed',PRT.model(m).model_name));
                end
            else
                beep;
                disp(sprintf('Model %s not properly specified! It will not be displayed',PRT.model(m).model_name));
            end
            
        end
        if ~nmi, error('There are no estimated/good models in this PRT!'); end
        
        handles.mi = mi;
        
        % Set model pulldown menu
        handles.mnames = model_name;
        set(handles.classmenu,'String',handles.mnames);
        
        % Get folds pulldown menu
        m             = get(handles.classmenu,'Value');
        handles.nfold = length(PRT.model(mi(m)).output.fold);
        folds{1}      = 'All folds / Average';
        for f = 1:handles.nfold
            folds{f+1} = num2str(f);
        end
        handles.folds = folds;
        set(handles.foldmenu,'String',handles.folds);
        
        % Set plots menu for first model
        if strcmp(PRT.model(mi(m)).input.type,'classification');
            if length(PRT.model(mi(m)).output.stats.c_acc) <= 2 ;
                plots = {'Histogram','Confusion Matrix','Predictions','ROC'};
            else
                plots = {'Histogram','Confusion Matrix','Predictions'};
            end
        else
            plots = {'Predictions (scatter)', 'Predictions (bar)', 'Predictions (line)'};
        end
        if isfield(PRT.model(mi(m)).input,'use_nested_cv')
            if PRT.model(mi(m)).input.use_nested_cv
                plots{length(plots)+1} = 'Influence of the hyper-parameter on performance';
            end
        end
        set(handles.plotmenu,'String',plots);
        
        % Initialize model button
        handles.model_button = 0;
        
        % Set the 'save permutations' weights' chackbox to 0
        handles.save_weights = 0;
        set(handles.save_perm_weights,'Value',0);
        %         set(handles.save_perm_weights,'Visible','off');
        %         set(handles.save_perm_weights,'Enable','off');
        
        % Clear axes
        cla(handles.axes5);
    end
end

% Choose default command line output for prt_ui_results
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_results wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_results_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in originbutton.
function originbutton_Callback(hObject, eventdata, handles)
% hObject    handle to originbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reset the crosshairs position
% -------------------------------------------------------------------------
if isfield(handles,'img')
    spm_orthviews('Reposition',[0 0 0]);
end

function mmedit_Callback(hObject, eventdata, handles)
% hObject    handle to mmedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of mmedit as text
%        str2double(get(hObject,'String')) returns contents of mmedit as a double

% Move crosshairs position in mm
% -------------------------------------------------------------------------
if isfield(handles,'img')
    mp    = handles.mmedit;
    posmm = get(mp,'String');
    pos = sscanf(posmm, '%g %g %g');
    if length(pos)~=3
        pos = spm_orthviews('Pos');
    end
    spm_orthviews('Reposition',pos);
end

% --- Executes during object creation, after setting all properties.
function mmedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mmedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function vxedit_Callback(hObject, eventdata, handles)
% hObject    handle to vxedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vxedit as text
%        str2double(get(hObject,'String')) returns contents of vxedit as a double

% Move crosshairs position in vx
% -------------------------------------------------------------------------
if isfield(handles,'img')
    mp    = handles.vxedit;
    posvx = get(mp,'String');
    pos   = sscanf(posvx, '%g %g %g');
    if length(pos)~=3
        pos = spm_orthviews('pos',1);
    end
    tmp = handles.vols{1}.mat;
    pos = tmp(1:3,:)*[pos ; 1];
    spm_orthviews('Reposition',pos);
end


% --- Executes during object creation, after setting all properties.
function vxedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vxedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in helpbutton.
function helpbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Help window for PRoNTo results has been launched.')
prt_ui_results_help;

% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close and clear figure
% -------------------------------------------------------------------------
close(handles.figure1);


function loadweight_Callback(hObject, eventdata, handles)
% hObject    handle to loadweight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadweight as text
%        str2double(get(hObject,'String')) returns contents of loadweight as a double


% --- Executes during object creation, after setting all properties.
function loadweight_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadweight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in weightbutton.
function weightbutton_Callback(hObject, eventdata, handles)
% hObject    handle to weightbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disp('Loading weights...>>')
% Select results (.img) for weight map
% -------------------------------------------------------------------------
if ~isfield(handles,'wmap') || ~handles.noloadw
    wmap            = spm_select(1,'image','Select weight map.');
    % Remove number in file name
    if strcmp(wmap(end-1),',')
        wmap = wmap(1:end-2);
    end
    V               = spm_vol(wmap);
    handles.vols{1} = V;
    handles.wmap    = wmap;
    
end

spm_orthviews('Reset');
if isfield(handles,'aimg')
    anatomicalbutton_Callback(hObject, eventdata, handles);
end

% Image dimensions
% -------------------------------------------------------------------------
fold          = get(handles.foldmenu,'Value')-1;
Vfolds        = handles.vols{1};
V             = Vfolds(1);
M             = V.mat;
DIM           = V.dim(1:3)';
xdim          = DIM(1); ydim  = DIM(2); zdim  = DIM(3);
if length(V.private.dat.dim) < 4
    fdim = 1;                       % Handle 3D images
else
    fdim = V.private.dat.dim(4);    % Handle 4D images
end
[xords,yords] = ndgrid(1:xdim,1:ydim);
xords         = xords(:)';  yords = yords(:)';
I             = 1:xdim*ydim;
zords_init    = ones(1,xdim*ydim);

% Get image values above zero for each fold and all folds
% -------------------------------------------------------------------------
xyz_above = [];
z_above   = [];
if fold == 0,
    fold_coord = fdim*ones(1,xdim*ydim);
    V = Vfolds(fdim);
else
    fold_coord = fold*ones(1,xdim*ydim);
    V = Vfolds(fold);
end

for z = 1:zdim,
    zords = z*zords_init;
    xyz   = [xords(I); yords(I); zords(I); fold_coord];
    zvals = spm_get_data(V,xyz);
    above = find(~isnan(zvals));
    if length(above)==length(zvals) %old version of weight computation
        above = find(zvals~=0);
    end
    if ~isempty(above)
        xyz_above = [xyz_above,xyz(:,above)];
        z_above   = [z_above,zvals(above)];
    end
end
XYZ   = xyz_above(1:3,:);
Z     = z_above;

% Set spm_orthviews properties
% -------------------------------------------------------------------------
rotate3d off
global st

handles.notinit = 1;
handles.img     = 1;

st.handles  = handles;
st.fig      = handles.figure1;
st.V        = V;
st.callback = 'prt_ui_results(''showpos'')';

% Display maps
% -------------------------------------------------------------------------
h  = spm_orthviews('Image', handles.wmap,[0.0619 0.0859 0.40 0.4196]);
handles.wimgh = h;
spm_orthviews('AddContext', h);
spm_orthviews('MaxBB');
spm_orthviews('AddBlobs', h, XYZ, Z, M);
cmap = get(gcf,'Colormap');
if size(cmap,1)~=128
    spm_figure('Colormap','jet');
end
spm_orthviews('Redraw');

% Show positions
% -------------------------------------------------------------------------
prt_ui_results('showpos');

disp('Done');

% Reset flag to load weights
handles.noloadw = 1;

% Show file name
% -------------------------------------------------------------------------
set(handles.loadweight,'String',handles.wmap);
guidata(hObject, handles);

function loadanatomical_Callback(hObject, eventdata, handles)
% hObject    handle to loadanatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of loadanatomical as text
%        str2double(get(hObject,'String')) returns contents of loadanatomical as a double


% --- Executes during object creation, after setting all properties.
function loadanatomical_CreateFcn(hObject, eventdata, handles)
% hObject    handle to loadanatomical (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in anatomicalbutton.
function anatomicalbutton_Callback(hObject, eventdata, handles)
% hObject    handle to anatomicalbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check if weight map and anatomical image exist and reset orthviews
% -------------------------------------------------------------------------
if ~isfield(handles,'wmap')
    spm_orthviews('Reset');
end
global st
st.fig = handles.figure1;
if ~isfield(handles,'aimg') || ~handles.noloadi
    img    = spm_select(1,'image','Select anatomical image.');
else
    img = handles.aimg;
end

% Show anatomical image
% -------------------------------------------------------------------------
rotate3d off
st.fig = handles.figure1;
handle = spm_orthviews('Image', img, [0.5295 0.0859 0.4196 0.4196]);
cmap   = get(gcf,'Colormap');
if size(cmap,1)~=128
    spm_figure('Colormap','gray')
end

handles.aimgh   = handle;
handles.aimg    = img;
handles.img     = 1;
handles.noloadi = 1;

% Show file name
% -------------------------------------------------------------------------
set(handles.loadanatomical,'String',handles.aimg);

guidata(hObject, handles);


% --- Executes on selection change in foldmenu.
function foldmenu_Callback(hObject, eventdata, handles)
% hObject    handle to foldmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns foldmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from foldmenu

% Change weight map
% -------------------------------------------------------------------------
if ~handles.model_button
    if isfield(handles,'vols')
        handles.noloadw = 1;
        weightbutton_Callback(hObject, eventdata, handles);
    end
end

% Change plot
% -------------------------------------------------------------------------
if isfield(handles,'plot')
    plotmenu_Callback(hObject, eventdata, handles);
end

% Change stats
% -------------------------------------------------------------------------
if isfield(handles,'stats')
    statsbutton_Callback(hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function foldmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in plotmenu.
function plotmenu_Callback(hObject, eventdata, handles)
% hObject    handle to plotmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns plotmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from plotmenu

plotm         = get(handles.plotmenu,'Value');
plotchosen    = num2str(plotm);
fold          = get(handles.foldmenu,'Value');
model         = get(handles.classmenu,'Value');
mi            = handles.mi;
model         = mi(model);
PRT           = handles.PRT;
handles.plot  = 1;
nms           = 7;
isyc1         = 0;
isyc2         = 0;
c1            = 0;
c2            = 0;
rotate3d off

if strcmp(PRT.model(model).input.type,'classification')
    %     mclass        = length(handles.PRT.model(model).output.stats.c_acc);
    %     if mclass ~= 2
    %         multiplot  = 4;
    %         plotchosen = num2str(multiplot(plotm));
    %     end
    
    % Check if the model has hyper-parameter optimization
    if (~isfield(handles.PRT.model(model).input,'use_nested_cv') || ...
            ~handles.PRT.model(model).input.use_nested_cv) && plotchosen == '4'
        plotchosen = '1';
    end
    
    % Plot
    % ---------------------------------------------------------------------
    switch plotchosen
        
        
        % Histograms
        % -----------------------------------------------------------------
        case '1'
            prt_plot_histograms(handles.PRT, model, fold, handles.axes5);
            
            % Confusion matrix
            % -------------------------------------------------------------
        case '2'
            prt_plot_confusion_matrix(handles.PRT, model, fold, handles.axes5);
            
            % Predictions
            % -------------------------------------------------------------
        case '3'
            prt_plot_prediction(handles.PRT, model, fold, nms, handles.axes5);
            
            % ROC / AUC
            % -------------------------------------------------------------
        case '4'
            prt_plot_ROC(handles.PRT, model, fold, handles.axes5);
            
            % TODO: Check if this does not cause problems when the
            % Influence of the hyper-parameter on performance was not used
        case '5'
            prt_plot_nested_cv(handles.PRT, model, fold, handles.axes5);
            
    end
    
else
    
    % Plot
    % ---------------------------------------------------------------------
    switch plotchosen
        case '1'
            prt_plot_prediction_reg_scatter(handles.PRT, model, handles.axes5);
            
        case '2'
            prt_plot_prediction_reg_bar(handles.PRT, model, handles.axes5);
            
        case '3'
            prt_plot_prediction_reg_line(handles.PRT, model, handles.axes5);
            
            % TODO: Check if this does not cause problems when the
            % nested CV was not used
        case '4'
            prt_plot_nested_cv(handles.PRT, model, fold, handles.axes5);
            
    end
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function plotmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in classmenu.
function classmenu_Callback(hObject, eventdata, handles)
% hObject    handle to classmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns classmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classmenu

% Hints: contents = cellstr(get(hObject,'String')) returns classmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classmenu

% Get folds
m  = get(handles.classmenu,'Value');
mi = handles.mi;
if length(handles.PRT.model(mi(m)).output)>1
    beep
    disp('Cannot display results per kernel in results window')
    return
end
handles.nfold = length(handles.PRT.model(mi(m)).output.fold);


folds{1}      = 'All folds / Average';
for f = 1:handles.nfold
    folds{f+1} = num2str(f);
end

% Set plots menu for first model
if strcmp(handles.PRT.model(mi(m)).input.type,'classification');
    if length(handles.PRT.model(mi(m)).output(1).stats.c_acc) <= 2 ;
        plots = {'Histogram','Confusion Matrix','Predictions','ROC'};
    else
        plots = {'Histogram', 'Confusion Matrix'};
    end
else
    plots = {'Predictions (scatter)', 'Predictions (bar)', 'Predictions (line)'};
end
if isfield(handles.PRT.model(mi(m)).input,'use_nested_cv')
    if handles.PRT.model(mi(m)).input.use_nested_cv
        plots{length(plots)+1} = 'Influence of the hyper-parameter on performance';
    end
end

set(handles.plotmenu,'String',plots);

% Set folds and call fold function to change plot/stats
handles.folds = folds;
set(handles.foldmenu,'String',handles.folds);
handles.model_button = 1;

foldmenu_Callback(hObject, eventdata, handles);

% Update stats if they are being shown
if isfield(handles, 'stats')
    statsbutton_Callback(hObject, eventdata, handles);
end

handles.model_button = 0;
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function classmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to classmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function repedit_Callback(hObject, eventdata, handles)
% hObject    handle to repedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of repedit as text
%        str2double(get(hObject,'String')) returns contents of repedit as a double


% --- Executes during object creation, after setting all properties.
function repedit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to repedit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in permutbutton.
function permutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to permutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

m    = get(handles.classmenu,'Value');
mi   = handles.mi;
reps = str2num(get(handles.repedit,'String'));
if  ~isempty(reps)
    if length(reps) ==1
        reps = round(reps);
        disp('Performing permutation test.........>>')
        prt_permutation(handles.PRT, reps, mi(m), handles.pathdir,...
            handles.save_weights);
        % Load new PRT.mat
        PRTmat = fullfile(handles.pathdir,'PRT.mat');
        load(PRTmat);
        perm     = PRT.model(mi(m)).output.stats.permutation;
        stats.perm      = perm;
        stats.show_perm = 1;
        if strcmp(handles.PRT.model(mi(m)).input.type,'classification');
            stats.type = 'class';
        else
            stats.type = 'reg';
        end
        handles.stats = stats;
        handles.PRT   = PRT;
        
        % Update GUI
        guidata(hObject, handles);
        
        % Call stats button
        statsbutton_Callback(hObject, eventdata, handles);
    else
        beep;
        disp('Please enter only one value for the number of repetitions!');
    end
else
    beep;
    disp('Repetitions should be a number!');
end

% --- Executes on button press in save_perm_weights.
function save_perm_weights_Callback(hObject, eventdata, handles)
% hObject    handle to save_perm_weights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of save_perm_weights
% Save the weights and predictions for each permutation if required
flag=get(handles.save_perm_weights,'Value');
if flag
    handles.save_weights=1;
else
    handles.save_weights=0;
end
guidata(hObject, handles);


% --- Executes on button press in statsbutton.
function statsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to statsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Reads model and fold
% -------------------------------------------------------------------------
fold  = get(handles.foldmenu,'Value');
mi    = handles.mi;
m     = get(handles.classmenu,'Value');
PRT   = handles.PRT;

% Check if stats exist
% -------------------------------------------------------------------------
if isfield(handles, 'stats')
    stats = handles.stats;
end

% Read stats
% -------------------------------------------------------------------------
if strcmp(PRT.model(mi(m)).input.type,'classification')
    if fold == 1
        macc  = PRT.model(mi(m)).output.stats.acc;  % overall acc
        mbacc = PRT.model(mi(m)).output.stats.b_acc;
        mcacc = PRT.model(mi(m)).output.stats.c_acc;
        mcpv  = PRT.model(mi(m)).output.stats.c_pv;
        if isfield(PRT.model(mi(m)).output.stats,'permutation') && ...
                ~isempty(PRT.model(mi(m)).output.stats.permutation)
            stats.show_perm=1;
            stats.perm.pvalue_b_acc=PRT.model(mi(m)).output.stats.permutation.pvalue_b_acc;
            stats.perm.pvalue_c_acc=PRT.model(mi(m)).output.stats.permutation.pvalue_c_acc;
        end
    else
        macc  = PRT.model(mi(m)).output.fold(fold-1).stats.acc;
        mbacc = PRT.model(mi(m)).output.fold(fold-1).stats.b_acc;
        mcacc = PRT.model(mi(m)).output.fold(fold-1).stats.c_acc;
        mcpv  = PRT.model(mi(m)).output.fold(fold-1).stats.c_pv;
    end
    
    stats.macc  = macc;
    stats.mbacc = mbacc;
    stats.mcacc = mcacc;
    stats.mcpv  = mcpv;
    stats.type  = 'class';
    
    prt_ui_stats(stats,handles.prtdir);
    
else
    if fold == 1
        corr  = PRT.model(mi(m)).output.stats.corr;  % overall correlation
        if isfield(PRT.model(mi(m)).output.stats,'r2')
            r2  = PRT.model(mi(m)).output.stats.r2;  % squared correlation
        end
        mse   = PRT.model(mi(m)).output.stats.mse;  % overall mse
        if isfield(PRT.model(mi(m)).output.stats,'nmse')
            nmse  = PRT.model(mi(m)).output.stats.nmse;  % normalised mse
        end
        if isfield(PRT.model(mi(m)).output.stats,'permutation') && ...
                ~isempty(PRT.model(mi(m)).output.stats.permutation)
            stats.show_perm=1;
            stats.perm.pval_corr = PRT.model(mi(m)).output.stats.permutation.pval_corr;
            stats.perm.pval_mse = PRT.model(mi(m)).output.stats.permutation.pval_mse;
            
            if isfield(PRT.model(mi(m)).output.stats.permutation, 'pval_r2'),...
                    stats.perm.pval_r2 = PRT.model(mi(m)).output.stats.permutation.pval_r2; end
            if isfield(PRT.model(mi(m)).output.stats.permutation, 'pval_nmse'),...
                    stats.perm.pval_nmse = PRT.model(mi(m)).output.stats.permutation.pval_nmse; end
            
        end
    else
        corr  = PRT.model(mi(m)).output.fold(fold-1).stats.corr;  % overall correlation
        if isfield(PRT.model(mi(m)).output.stats,'r2')
            r2  = PRT.model(mi(m)).output.fold(fold-1).stats.r2;  % overall correlation
        end
        mse   = PRT.model(mi(m)).output.fold(fold-1).stats.mse;  % overall mse
        if isfield(PRT.model(mi(m)).output.stats,'nmse')
            nmse  = PRT.model(mi(m)).output.fold(fold-1).stats.nmse;  % normalised mse
        end
    end
    
    stats.corr = corr;
    if isfield(PRT.model(mi(m)).output.stats,'r2'), stats.r2 = r2; end
    stats.mse  = mse;
    if isfield(PRT.model(mi(m)).output.stats,'nmse'), stats.nmse = nmse; end
    stats.type = 'reg';
    
    prt_ui_stats(stats,handles.prtdir);
    
end

stats.show_perm = 0;
handles.stats = stats;
guidata(hObject, handles);

% Show crosshairs position
% -------------------------------------------------------------------------
function showpos()

global st

mp13 = st.handles.mmedit;
mp14 = st.handles.vxedit;
tx20 = st.handles.posintensitytext;

set(mp13,'String',sprintf('%.1f %.1f %.1f',spm_orthviews('Pos')));
pos = spm_orthviews('Pos',1);
set(mp14,'String',sprintf('%.1f %.1f %.1f',pos));
set(tx20,'String',sprintf('%g',spm_sample_vol(st.V,pos(1),pos(2),pos(3),st.hld)));

cmap = get(gcf,'Colormap');
if size(cmap,1)~=128
    spm_figure('Colormap','gray-jet');
end


% --- Executes on button press in resetbutton.
function resetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to resetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
spm_orthviews('Reset');
if isfield(handles, 'wmap'), handles = rmfield(handles, 'wmap'); end
if isfield(handles, 'aimg'), handles = rmfield(handles,'aimg'); end
handles.noloadw = 0;
guidata(hObject, handles);

% Save menu
% -------------------------------------------------------------------------
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
    'Save figure as','.png');
[a,b,c]=fileparts(filename);
ext=['-d',c(2:end)];

% Set the color of the different backgrounds and figure parameters to white
cf=get(handles.figure1,'Color');
set(handles.figure1,'Color',[1,1,1])
aa=get(handles.figure1,'children');
xc=[];
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        try
            xc=[xc;get(aa(i),'BackgroundColor')];
            set(aa(i),'BackgroundColor',[1 1 1])
        end
        bb=get(aa(i),'children');
        if ~isempty(bb)
            for j=1:length(bb)
                try
                    xc=[xc;get(bb(j),'BackgroundColor')];
                    set(bb(j),'BackgroundColor',[1 1 1])
                end
                if strcmpi(get(bb(j),'type'),'uipanel')
                    cc=get(bb(j),'children');
                    if ~isempty(cc)
                        for k=1:length(cc)
                            try
                                xc=[xc;get(cc(k),'BackgroundColor')];
                                set(cc(k),'BackgroundColor',[1 1 1])
                            end
                            if strcmpi(get(cc(k),'type'),'uipanel')
                                dd=get(cc(k),'children');
                                if ~isempty(dd)
                                    for l=1:length(dd)
                                        try
                                            xc=[xc;get(dd(l),'BackgroundColor')];
                                            set(dd(l),'BackgroundColor',[1 1 1])
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            xc=[xc;get(aa(i),'BackgroundColor')];
            set(aa(i),'BackgroundColor',[1 1 1])
        end
    end
end

if ~strcmp(ext, '-dfig')
    print(handles.figure1,ext,[pathname,filesep,b],'-r500')
else
    saveas(handles.figure1,[pathname,filesep,b],'fig')
end


% Set the color of the different backgrounds and figure parameters to white
set(handles.figure1,'Color',cf)
scount=1;
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        try
            set(aa(i),'BackgroundColor',xc(scount,:))
            scount=scount+1;
        end
        bb=get(aa(i),'children');
        if ~isempty(bb)
            for j=1:length(bb)
                try
                    set(bb(j),'BackgroundColor',xc(scount,:))
                    scount=scount+1;
                end
                if strcmpi(get(bb(j),'type'),'uipanel')
                    cc=get(bb(j),'children');
                    if ~isempty(cc)
                        for k=1:length(cc)
                            try
                                set(cc(k),'BackgroundColor',xc(scount,:))
                                scount=scount+1;
                            end
                            if strcmpi(get(cc(k),'type'),'uipanel')
                                dd=get(cc(k),'children');
                                if ~isempty(dd)
                                    for l=1:length(dd)
                                        try
                                            set(dd(l),'BackgroundColor',xc(scount,:))
                                            scount=scount+1;
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    elseif ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            set(aa(i),'BackgroundColor',xc(scount,:))
            scount=scount+1;
        end
    end
end

cd(wd)
