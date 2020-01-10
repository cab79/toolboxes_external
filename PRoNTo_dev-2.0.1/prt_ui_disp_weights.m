function varargout = prt_ui_disp_weights(varargin)
% PRT_UI_DISP_WEIGHTS MATLAB code for prt_ui_disp_weights.fig
%
% PRT_UI_DISP_WEIGHTS, by itself, creates a new PRT_UI_DISP_WEIGHTS or raises the
% existing singleton*.
%
% H = PRT_UI_DISP_WEIGHTS returns the handle to a new PRT_UI_DISP_WEIGHTS or the
% handle to the existing singleton*.
%
% PRT_UI_DISP_WEIGHTS('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_DISP_WEIGHTS.M with the given input arguments.
%
% PRT_UI_DISP_WEIGHTS('Property','Value',...) creates a new PRT_UI_DISP_WEIGHTS or
% raises the existing singleton*.  Starting from the left, property value
% pairs are applied to the GUI before prt_ui_disp_weights_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_disp_weights_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa and J. Schrouff
% $Id: prt_ui_disp_weights.m 784 2013-08-22 16:43:32Z monteiro $

% Edit the above text to modify the response to help prt_ui_disp_weights

% Last Modified by GUIDE v2.5 09-Mar-2015 16:15:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prt_ui_disp_weights_OpeningFcn, ...
    'gui_OutputFcn',  @prt_ui_disp_weights_OutputFcn, ...
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


% --- Executes just before prt_ui_disp_weights is made visible.
function prt_ui_disp_weights_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_disp_weights (see VARARGIN)

%if window already exists, just put it as the current figure
Tag='WeightResults';
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
    set(handles.figure1,'Name','PRoNTo :: Weights')
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
        if isempty(PRT)
            error('prt_ui_disp_weights:NoPRT','No PRT file selected')        
        end

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
        nmroi = 0;
        for m = 1:nmodels
            if isfield(PRT.model(m),'input') && ~isempty(PRT.model(m).input)
                if isfield(PRT.model(m),'output') && ~isempty(PRT.model(m).output)
                    if isfield(PRT.model(m).output,'weight_img') && ~isempty(PRT.model(m).output.weight_img)
                        nmi = nmi +1;
                        model_name{nmi} = PRT.model(m).model_name;
                        mi = [mi, m];
                        if isfield(PRT.model(m).output,'weight_ROI') && ...
                                ~isempty(PRT.model(m).output.weight_ROI)
                            nmroi = nmroi+1;
                        end
                    else
                        disp(sprintf('Weights not computed for model %s ! It will not be displayed',PRT.model(m).model_name));
                    end
                else
                    beep;
                    disp(sprintf('Weights not computed for model %s ! It will not be displayed',PRT.model(m).model_name));
                end
            else
                beep;
                disp(sprintf('Model %s not properly specified! It will not be displayed',PRT.model(m).model_name));
            end
            
        end
        %         if ~nmi, error('There are no estimated/good models in this
        %         PRT!'); end
        
        if nmi
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
            set(handles.foldmenu,'Value',1);
            
            list = handles.PRT.model(mi(m)).output.weight_img;
            if ~isempty(list)
                set(handles.imagemenu,'String',list)
                set(handles.imagemenu,'Value',1)
            else
                set(handles.imagemenu,'Enable','off')
            end
            
            
        end
        % Initialize model button
        handles.model_button = 0;
        
        % Deal with ROI table and bar graph: 1st turned off        
        set(handles.axes1, 'visible','off')
        set(handles.ROItable,'visible','off')
        handles.selectedcell = [];
        handles.labels=cell(nmodels,1);
        set(handles.butt_load_labels,'visible','off');
        set(handles.saveweight,'visible','off');
        set(handles.modtable,'Visible','off')
        set(handles.modtable,'Enable','off')
        % Clear axes
        cla(handles.axes1);
%         classmenu_Callback(hObject, eventdata, handles)
    end
end
set(handles.disp_voxels,'Value',1)
set(handles.disp_voxels,'Enable','off')
set(handles.disp_regions,'Value',0)
set(handles.disp_regions,'Enable','off')
uistack(handles.weightspanel,'bottom')
% Choose default command line output for prt_ui_disp_weights
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_disp_weights wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_disp_weights_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
% delete(handles.figure1)
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
    child = get(handles.weightspanel,'Children');
    count = [];
    for i = 1:length(child) % Resize colorbar if blobs
        if strcmpi(get(child(i),'Type'),'axes')
            count = [count, i];
        end
    end
    if length(count)==4
        pos = get(child(count(1)),'Position');
        set(child(count(1)),'Position',[pos(1)*1.3,pos(2),pos(3),pos(4)*0.9])
    end
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
    child = get(handles.weightspanel,'Children');
    count = [];
    for i = 1:length(child) % Resize colorbar if blobs
        if strcmpi(get(child(i),'Type'),'axes')
            count = [count, i];
        end
    end
    if length(count)==4
        pos = get(child(count(1)),'Position');
        set(child(count(1)),'Position',[pos(1)*1.3,pos(2),pos(3),pos(4)*0.9])
    end
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
    child = get(handles.weightspanel,'Children');
    count = [];
    for i = 1:length(child) % Resize colorbar if blobs
        if strcmpi(get(child(i),'Type'),'axes')
            count = [count, i];
        end
    end
    if length(count)==4
        pos = get(child(count(1)),'Position');
        set(child(count(1)),'Position',[pos(1)*1.3,pos(2),pos(3),pos(4)*0.9])
    end
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


% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Close and clear figure
% -------------------------------------------------------------------------
Tag='WeightResults';
F = findall(allchild(0),'Flat','Tag',Tag);
delete(F);
clear handles
clear F


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
function handles = weightbutton_Callback(hObject, eventdata, handles)
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
gcfchil = get(gcf,'Children');
for i=1:length(gcfchil)
    if strcmpi(get(gcfchil(i),'Type'),'axes')
        delete(gcfchil(i)) % Delete any left-over colorbar
    end
end
gcfchil = get(handles.weightspanel,'Children');
for i=1:length(gcfchil)
    if strcmpi(get(gcfchil(i),'Type'),'axes')
        delete(gcfchil(i)) % Delete any left-over colorbar
    end
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

if ~isempty(handles.selectedcell)
    matroi = NaN * ones(xdim,ydim,zdim);
    idroi = handles.idfeat_roi{handles.sort_roi(handles.selectedcell(1))};
    matroi(idroi) = 1;
    handles.roimatdisp = matroi;
else
    matroi = ones(xdim,ydim,zdim);
    handles.roimatdisp = matroi;
end

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
    if isfield(handles,'roimatdisp') && ~isempty(handles.roimatdisp)
        indmask = sub2ind([xdim,ydim,zdim],xords(I)', yords(I)', zords(I)');
        mroi = handles.roimatdisp(indmask);
        zvals = zvals.*mroi';
    end
    above = find(~isnan(zvals));
    if length(above)==length(zvals) %old version of weight computation
        above = find(zvals~=0);
    end
    if ~isempty(above)
        xyz_above = [xyz_above,xyz(:,above)];
        z_above   = [z_above,zvals(above)];
    end
end
if ~isempty(xyz_above)
    XYZ   = xyz_above(1:3,:);
    Z     = z_above;
    
    %compute center of gravity to reposition crosshairs
    xm = round(median(XYZ(1,:)));
    ym = round(median(XYZ(2,:)));
    zm = round(median(XYZ(3,:)));
else
    disp('Only null contributions in this image')
end

% xmin = min(XYZ(1,:));
% xmax = max(XYZ(1,:));
% xfov =xmax-xmin;
% Set spm_orthviews properties
% -------------------------------------------------------------------------
rotate3d off
try
    clear st
end
global st

handles.notinit = 1;
handles.img     = 1;

st.handles  = handles;
st.fig      = handles.figure1;
st.V        = V;
st.callback = 'prt_ui_disp_weights(''showpos'')';

% Display maps
% -------------------------------------------------------------------------
[BB vx] = spm_get_bbox(handles.wmap);
xax = BB(1,1):abs(vx(1)):BB(2,1);
yax = BB(1,2):abs(vx(2)):BB(2,2);
zax = BB(1,3):abs(vx(3)):BB(2,3);

h  = spm_orthviews('Image', handles.wmap,[0.03 0.01 0.95 1.08]);
% h  = spm_orthviews('Image', handles.wmap,[0.03 0.01 0.95
% 1.08],handles.weightspanel); %to change when spm_orthviews will be
% updated, Matlab 2014b issue
handles.wimgh = h;
spm_orthviews('AddContext', h);
spm_orthviews('MaxBB');
if ~isempty(xyz_above)
    spm_orthviews('AddBlobs', h, XYZ, Z, M);
    cmap = get(gcf,'Colormap');
    if size(cmap,1)~=128
        spm_figure('Colormap','jet');
    end
    spm_orthviews('Reposition',[sign(vx(1))*xax(xm),sign(vx(2))*yax(ym),sign(vx(3))*zax(zm)])
    % spm_orthviews('Zoom',(xfov*abs(vx(1))))
    spm_orthviews('Redraw');
end

[fpw,faw]=fileparts(handles.wmap); %Sometimes extensions can pose problems
for i=1:length(st.vols) % Image was added to variable st
    if ~isempty(st.vols{i})
        [fp,fa] = fileparts(st.vols{i}.fname);
        if strcmpi([fp,fa],[fpw,faw])
            idxw = i;
            set(st.vols{idxw}.ax{1}.ax,'parent',handles.weightspanel);
            set(st.vols{idxw}.ax{2}.ax,'parent',handles.weightspanel);
            set(st.vols{idxw}.ax{3}.ax,'parent',handles.weightspanel);
            if ~isempty(xyz_above)
                set(st.vols{idxw}.blobs{1}.cbar,'parent',handles.weightspanel);
                cbp = get(st.vols{idxw}.blobs{1}.cbar,'Position'); % Colorbar position
                set(st.vols{idxw}.blobs{1}.cbar,'Position',...
                    [cbp(1)*1.3,cbp(2),cbp(3),cbp(4)*0.9]);
                handles.posaxe4 = get(st.vols{idxw}.blobs{1}.cbar,'Position');
            end
            handles.posaxe1 = get(st.vols{idxw}.ax{1}.ax,'Position');
            handles.posaxe2 = get(st.vols{idxw}.ax{2}.ax,'Position');
            handles.posaxe3 = get(st.vols{idxw}.ax{3}.ax,'Position');
        end
    end
end

if isfield(handles,'aimg')
    anatomicalbutton_Callback(hObject, eventdata, handles);
end

% Show positions
% -------------------------------------------------------------------------
prt_ui_disp_weights('showpos');

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
% gcfchil = get(handles.anatomicalpanel,'Children');
% for i=1:length(gcfchil)
%     if strcmpi(get(gcfchil(i),'Type'),'axes')
%         delete(gcfchil(i)) % Delete any left-over colorbar
%     end
% end
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
handle = spm_orthviews('Image', img,[0.03 0.01 0.95 1.08]);
% handle = spm_orthviews('Image', img,[0.03 0.01 0.95 1.08],handles.anatomicalpanel);
spm_orthviews('AddContext',handle);
[pimg,aimg]=fileparts(img); % Sometimes extension can pose problems
if isfield(handles,'wmap') 
    [fpl,fal] = fileparts(handles.wmap);
end
for i=1:length(st.vols) % Image was added to variable st
    if ~isempty(st.vols{i}) % Look for weight image
        [fp,fa] = fileparts(st.vols{i}.fname);
        if isfield(handles,'wmap') && strcmpi([fp,fa],[fpl,fal])
            idxw = i;
            set(st.vols{idxw}.ax{1}.ax,'Position',handles.posaxe1);
            set(st.vols{idxw}.ax{2}.ax,'Position',handles.posaxe2);
            set(st.vols{idxw}.ax{3}.ax,'Position',handles.posaxe3);
            set(st.vols{idxw}.blobs{1}.cbar,'Position',handles.posaxe4);
        elseif strcmpi([fp,fa],[pimg,aimg]) % Look for anatomical image
            idxa = i;
            set(st.vols{idxa}.ax{1}.ax,'parent',handles.anatomicalpanel);
            set(st.vols{idxa}.ax{2}.ax,'parent',handles.anatomicalpanel);
            set(st.vols{idxa}.ax{3}.ax,'parent',handles.anatomicalpanel);
        end
    end
end

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


% Update table and bar graph
%--------------------------------------------------------------------------
m  = get(handles.classmenu,'Value');
if m==0
    m=1;
end
mi = handles.mi;
ffi = get(handles.foldmenu,'Value')-1;
if ffi == -1 % for Mac issues with popup menus
    ffi = 0;
end
if ffi==0 % for average
    ffi=length(get(handles.foldmenu,'String'));
end
if isfield(handles.PRT.model(mi(m)).output,'weight_ROI') &&... % chosen model has ROI values
         ~isempty(handles.PRT.model(mi(m)).output.weight_ROI) &&...
       ~isempty(handles.PRT.model(mi(m)).output.weight_ROI{1})
   dat = handles.dattable;
   weights = handles.PRT.model(mi(m)).output.weight_ROI{handles.class}(:,ffi)*100;
   dat(:,2) = num2cell(weights);
   [vald,idwroi] = sort(weights,'descend');
   dat = dat(idwroi,:);
   handles.sort_roi= idwroi;
   set(handles.ROItable,'Data',dat);
   set(handles.ROItable,'visible','on');
   set(handles.butt_load_labels,'visible','on');
   
   %Bar graph to show decrease in ROI weights
   set(handles.axes1,'visible','on')
   bar(handles.axes1,weights(idwroi));
   set(get(handles.axes1,'XLabel'),'FontWeight','demi')
   set(get(handles.axes1,'XLabel'),'String','ROI index')
   set(get(handles.axes1,'YLabel'),'String','ROI weight')
   set(get(handles.axes1,'YLabel'),'FontWeight','demi')    
   
   if ~isempty(handles.datmod)
       datmod = handles.datmod;
       wemod = zeros(length(handles.nmods),handles.nfold+1);
        for i=1:length(handles.nmods) % get all the modalities
            wemod(i,:) = [handles.PRT.model(mi(m)).output.weight_MOD{i}];
        end
        datmod(:,2) = num2cell(wemod(:,ffi)*100);
        [vald,idwmod] = sort(wemod(:,ffi),'descend');
        datmod= datmod(idwmod,:);
        set(handles.modtable,'Data',datmod)
        set(handles.modtable,'Visible','on')
   end
end

% Change weight map
% -------------------------------------------------------------------------
if ~handles.model_button
    if isfield(handles,'vols')
        handles.noloadw = 1;
        handles.selectedcell = [];
        weightbutton_Callback(hObject, eventdata, handles);
    end
end

guidata(hObject, handles);




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


% --- Executes on selection change in classmenu.
function classmenu_Callback(hObject, eventdata, handles)
% hObject    handle to classmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns classmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classmenu

% Hints: contents = cellstr(get(hObject,'String')) returns classmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from classmenu
m  = get(handles.classmenu,'Value');
if m==0
    m=1;
end
mi = handles.mi;

% Set the list of weight images generated for this model
list = handles.PRT.model(mi(m)).output.weight_img;
if ~isempty(list)
    set(handles.imagemenu,'String',list)
else
    set(handles.imagemenu,'Enable','off')
end
set(handles.imagemenu,'Value',1);
handles.class = 1;

% Get the modalities in the model if multiple kernels and multiple
% modalities
in.fs_name = handles.PRT.model(mi(m)).input.fs(1).fs_name;
fid = prt_init_fs(handles.PRT,in);
handles.fid = fid;
if handles.PRT.fs(fid).multkernel
    nmod = length(handles.PRT.fs(fid).modality);
    mods = cell(nmod,1);
    for i=1:nmod
        mods{i} = handles.PRT.fs(fid).modality(i).mod_name;
    end
    if ~isempty(strfind(handles.PRT.model(mi(m)).input.machine.function,'MKL'))
        handles.summed = 0;
    else
        handles.summed = 1; % summed modalities
    end
elseif ~isempty(strfind(handles.PRT.model(mi(m)).input.machine.function,'MKL')) && ...
        ~handles.PRT.fs(fid).multkernel
    mods{1} = handles.PRT.fs(fid).modality(1).mod_name;
    handles.summed = 0; 
else
    mods{1} = handles.PRT.fs(fid).modality(1).mod_name;
    handles.summed = 1; % concatenated or one modality
end
handles.nmods = mods;
guidata(hObject, handles);
imagemenu_Callback(hObject,eventdata,handles)
handles = guidata(hObject);
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

% --- Executes on selection change in imagemenu.
function imagemenu_Callback(hObject, eventdata, handles)
% hObject    handle to imagemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns imagemenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from imagemenu
% Get models (classmenu), classes (imagemenu) and folds (foldmenu)
m  = get(handles.classmenu,'Value');
if m==0
    m=1;
end
mi = handles.mi;
handles.nfold = length(handles.PRT.model(mi(m)).output.fold);


folds{1}      = 'All folds / Average';
for f = 1:handles.nfold
    folds{f+1} = num2str(f);
end

% Set folds and call fold function to plot the weights for that model
handles.folds = folds;
set(handles.foldmenu,'String',handles.folds);

disp_vox = get(handles.disp_voxels,'Value');
set(handles.disp_voxels,'Enable','on')

%get the fold to display
ffi = get(handles.foldmenu,'Value')-1;
if ffi == -1 % for Mac issues with popup menus
    ffi = 0;
end
if ffi==0 % for average
    ffi=length(get(handles.foldmenu,'String'));
end

% get the image to display
handles.class = get(handles.imagemenu,'Value');
if handles.class ==0
    handles.class = 1;
end

% Compare the name of the image with the name of the modalities to obtain
% the index of the modality
nim = get(handles.imagemenu,'String');
for i = 1:length(handles.nmods)
    if ~isempty(strfind(nim(handles.class),char(handles.nmods{i})))
        mids = i;
    end
end


% Initialize the table
% -------------------------------------------------------------------------

flagmodMKL = 0; % One weight per modality?

% chosen model has ROI values or weights per modality : create the table
if isfield(handles.PRT.model(mi(m)).output,'weight_ROI') &&... 
        ~isempty(handles.PRT.model(mi(m)).output.weight_ROI) &&...
       ~isempty(handles.PRT.model(mi(m)).output.weight_ROI{1})
   
   in = struct();
   in.fs_name = handles.PRT.model(mi(m)).input.fs(1).fs_name;
   fid = prt_init_fs(handles.PRT,in);
   % For each modality if multiple modalities in multiple kernel settings
   if handles.PRT.fs(fid).multkernel 
       stm = length(handles.nmods); %if multiple modalities 
   else
       stm =1; % if modalities concatenated or only one
   end
   if ~iscell(handles.PRT.fs(fid).atlas_name)
       handles.PRT.fs(fid).atlas_name = {handles.PRT.fs(fid).atlas_name};
   end
   % Get the number of ROIs and their labels for each modality
   handles.num_roi = cell(stm,1);
   
   for i = 1:stm
       if isfield(handles.PRT.fs(fid).modality(i),'num_ROI') && ~handles.summed
            num_roi = handles.PRT.fs(fid).modality(i).num_ROI; % MKL on regions
            atl_name = handles.PRT.fs(fid).atlas_name{i};
       else
           num_roi = 1:size(handles.PRT.model(mi(m)).output.weight_ROI{handles.class},1);
           if isfield(handles.PRT.model(mi(m)).output,'weight_atlas') && ...
                   ~isempty(handles.PRT.model(mi(m)).output.weight_atlas)
                atl_name = handles.PRT.model(mi(m)).output.weight_atlas{i} ; %summarized weights
           else
               atl_name = []; %MKL on modalities
           end
       end
       handles.num_roi{i} = num_roi;
       
       % Get the labels if they are stored in a .mat along the atlas file for
       % MKL on regions and summarisation
       if ~isempty(atl_name)
           [a,b]=fileparts(atl_name);
           try
               load(fullfile(a,['Labels_',b,'.mat']))
               try
                   handles.labels{mi(m)}{i}=ROI_names;
               catch
                   disp('No variable ROI_names found, generic names used')
               end
           catch
                 disp('No file containing the names of the ROIs found, generic names used')
           end
           if ~isfield(handles,'labels') || ...
                   isempty(handles.labels{mi(m)}) || ...
                   i>length(handles.labels{mi(m)}) || ...
                   isempty(handles.labels{mi(m)}{i})
               label=cell(length(num_roi),1);
               for j=1:length(num_roi)
                   label{j} = ['ROI_',num2str(num_roi(j))];
               end
           else
               if isfield(handles.PRT.fs(fid),'multkernelROI') && ...
                       handles.PRT.fs(fid).multkernelROI && ...
                       ~handles.summed
                    if length(handles.PRT.fs(fid).modality(i).igood_kerns)==length(num_roi)
                        label = handles.labels{mi(m)}{i}(handles.PRT.fs(fid).modality(i).igood_kerns); %take 0 kernels out
                        handles.num_roi{i} = handles.PRT.fs(fid).modality(i).igood_kerns;
                    end
               else
                   try
                       label = handles.labels{mi(m)}{i}(num_roi);
                   catch
                       if length(num_roi)==length(handles.labels{mi(m)}{i})  % match between labels and number of regions
                           label = handles.labels{mi(m)}{i}(num_roi);
                           handles.num_roi{i} = num_roi;
                       elseif length(num_roi)==length(handles.labels{mi(m)}{i})+1    %missing the 'others' region
                           handles.labels{mi(m)}{i} = [{'others'};handles.labels{mi(m)}{i}];
                           label = handles.labels{mi(m)}{i}(num_roi);
                           handles.num_roi{i} = num_roi;
                       else
                           warning('prt_ui_disp_weights:LabelsDoNotCorrespondtoROI',... % could not use file
                               'Number of labels in .mat does not correspond to number of ROIs');
                           label=cell(length(num_roi(:,i)),1);
                           for j=1:length(handles.num_roi{i})
                               label{j} = ['ROI_',num2str(num_roi(j,i))];
                           end
                       end
                   end

               end
           end 
           lc = {'ROI label','ROI weight (%)'};
           xlabel = 'Index in table';
           ylabel = 'Weight';
           set(handles.butt_load_labels,'Visible','on')
           set(handles.saveweight,'Visible','on')
       else %get names of the modalities for MKL on modalities
           label = cell(length(handles.PRT.fs(fid).modality),1);
           for j = 1:length(handles.PRT.fs(fid).modality)
               label{j} = char(handles.PRT.fs(fid).modality(j).mod_name);
           end
           lc = {'Modality','Mod. weight (%)'};
           xlabel = 'Index in table';
           ylabel = 'Weight';
           set(handles.butt_load_labels,'Visible','off')
           set(handles.saveweight,'Visible','off')
       end
   end
   dat(:,1) = label;
   weights = handles.PRT.model(mi(m)).output.weight_ROI{handles.class}(:,ffi)*100;
   dat(:,2) = num2cell(weights);   
   
   % Fill the table with ROI size if ROIs
   if ~isempty(strfind(handles.PRT.model(mi(m)).input.machine.function,'MKL')) && ...
           handles.PRT.fs(fid).multkernelROI  %Multiple kernel learning on ROIs
       if isfield(handles.PRT.fs(fid).modality,'idfeat_img')&& ... % Get the indexes of each ROI in the image
               ~isempty(handles.PRT.fs(fid).modality(mids).idfeat_img)
           lc = [lc,{'ROI size (vox)'}];
           for i=1:length(handles.PRT.fas)
               if strcmpi(handles.PRT.fs(fid).modality(mids).mod_name,...
                       handles.PRT.fas(i).mod_name)
                   mid = i;
               end
           end
           idfeat = handles.PRT.fas(mid).idfeat_img;
           if isempty(handles.PRT.fs(fid).modality(mids).idfeat_fas) % get the 2nd level masking
               idfeat_fas = 1:length(idfeat);
           else
               idfeat_fas = handles.PRT.fs(fid).modality(mids).idfeat_fas;
           end
           handles.idfeat_roi = cell(length(handles.PRT.fs(fid).modality(mids).idfeat_img),1);
           for i = 1:length(handles.PRT.fs(fid).modality(mids).idfeat_img)
               dat(i,3) = {length(handles.PRT.fs(fid).modality(mids).idfeat_img{i})};
               handles.idfeat_roi{i} = idfeat(idfeat_fas(handles.PRT.fs(fid).modality(mids).idfeat_img{i}));
           end
       end
       set(handles.disp_regions,'Enable','on')
       flagmodMKL = 0;
   elseif isfield(handles.PRT.model(mi(m)).output,'weight_idfeatroi') &&... % Summarizing the weights for ROI
           ~isempty(handles.PRT.model(mi(m)).output.weight_idfeatroi) % Get the indexes of each ROI in the image
       lc = [lc,{'ROI size (vox)'}];
       handles.idfeat_roi = cell(length(handles.PRT.model(mi(m)).output.weight_idfeatroi{handles.class}),1);
       for i = 1:length(handles.PRT.model(mi(m)).output.weight_idfeatroi{handles.class})
           dat(i,3) = {length(handles.PRT.model(mi(m)).output.weight_idfeatroi{handles.class}{i})};
           handles.idfeat_roi{i} = handles.PRT.model(mi(m)).output.weight_idfeatroi{handles.class}{i};
       end
       set(handles.disp_regions,'Enable','on')
       flagmodMKL = 0;
   else % if MKL on modalities, do not fill third column and do not display weights per region
       set(handles.disp_regions,'Enable','off')
       flagmodMKL = 1;
   end 
   
   % Compute Expected Ranking for model
   w_all = handles.PRT.model(mi(m)).output.weight_ROI{handles.class}(:,1:handles.nfold)*100;
   %take 0 columns out of the computation
   mw = min(w_all);
   mxw = max(w_all);
   id0 = find(mw==mxw);
   id0 = id0(mw(id0)==0);
   id = 1:size(w_all,2);
   idtr = setdiff(id,id0);
   if isempty(idtr)   % if a modality has always 0 weights
       erwn=NaN*ones(length(num_roi),1);
   else
       w_all = w_all(:,idtr);
       %deal with NaNs
       [d1,d2]=sort(w_all,1,'descend');
       isn=find(isnan(w_all(:,1)));
       d3=1:length(isn);
       d4=length(isn)+1:size(d1,1);
       ihn=[d2(d4,:);d2(d3,:)];
       [d1,dwn]=sort(ihn);
       erwn=zeros(length(num_roi),1);
       for i=1:length(num_roi)
           for j=1:length(num_roi)
               tmp=length(find(dwn(i,1:end-1)==j));
               erwn(i)=erwn(i)+j*tmp;
           end
       end
       erwn=erwn/length(idtr);
   end
   
   dat = [dat, num2cell(erwn)]; 
   lc = [lc,{'Exp. Ranking'}];
   handles.dattable = dat;
   [vald,idwroi] = sort(weights,'descend');
   handles.sort_roi= idwroi;
else
    % Cut window
    set(handles.disp_regions,'Enable','off')
end

if ~iscell(handles.PRT.model(mi(m)).output.weight_img)
    handles.PRT.model(mi(m)).output.weight_img={handles.PRT.model(mi(m)).output.weight_img};
end

% Select the image to display and check whether its ROI version exists
if disp_vox
    fntl = handles.PRT.model(mi(m)).output.weight_img{handles.class};
    set(handles.disp_voxels,'Value',1)
    set(handles.disp_regions,'Value',0)
elseif ~disp_vox && isfield(handles.PRT.model(mi(m)).output,'weight_ROI') &&... % chosen model has ROI values
         ~isempty(handles.PRT.model(mi(m)).output.weight_ROI) &&...
       ~isempty(handles.PRT.model(mi(m)).output.weight_ROI{1})
    fntl = ['ROI_',handles.PRT.model(mi(m)).output.weight_img{handles.class}];%get only first image for now
    set(handles.disp_voxels,'Value',0)
    set(handles.disp_regions,'Value',1)
else
    beep
    disp('Weights per region selected, but no ROI_ image found')
    disp('Displaying image of weights per voxels')
    fntl = handles.PRT.model(mi(m)).output.weight_img{handles.class};
    set(handles.disp_voxels,'Value',1)
    set(handles.disp_regions,'Value',0)
end

% Call weightbutton to display the chosen image
if exist([handles.prtdir,filesep,fntl,'.img'],'file') || ...
        exist([handles.prtdir,filesep,fntl],'file')
    [a,b] = fileparts(fntl);
    handles.wmap = [handles.prtdir,filesep,b,'.img'];
    V               = spm_vol(handles.wmap);
    handles.vols{1} = V;
    handles.noloadw = 1;
    handles.model_button = 0;
    handles.selectedcell = [];
    % Update handles structure
    guidata(hObject, handles);
    handles=weightbutton_Callback(hObject, eventdata, handles);
end

% % Compute homogeneity of ROIs if MKL on ROIs or summarizing
% if ~flagmodMKL && isfield(handles.PRT.model(mi(m)).output,'weight_ROI') &&... % chosen model has ROI values
%        ~isempty(handles.PRT.model(mi(m)).output.weight_ROI) &&...
%         isfield(handles,'wmap') && ~isempty(handles.wmap)
%     
%     % Get homogeneity of each ROI in terms of sign from the image
%     dat = handles.dattable;
%     [vald,idwroi] = sort([dat{:,2}],'descend');
%     hom = zeros(length(vald),handles.nfold+1);
%     for f = 1:handles.nfold+1
%         imtl = V(f);
%         vv = spm_read_vols(imtl);
%         for i = 1:length(vald) %for each ROI
%             val = vv(handles.idfeat_roi{i});
%             hom(i,f) = length(find(val>0))/length(val);
%         end
%         if length(unique(hom(:,f)))==1 % probably all weights to zero
%             hom(:,f) = NaN;
%         end
%     end
%    handles.hom_roi = hom;
%    dat =[dat, num2cell(hom(:,end)*100)];
%    lc = [lc,{'# Pos. (%)'}];
% end

% Fill modality table in the case of multiple modalities and multiple ROIs
% (whether summarized or MKL)
if isfield(handles.PRT.model(mi(m)).output,'weight_MOD') &&...
        ~isempty(handles.PRT.model(mi(m)).output.weight_MOD{handles.class}) && ...
        ~flagmodMKL
    datmod = cell(length(handles.nmods),3);
    datmod(:,1) = [handles.nmods{:}];
    wemod = zeros(length(handles.nmods),handles.nfold+1);
    for i=1:length(handles.nmods) % get all the modalities
        wemod(i,:) = [handles.PRT.model(mi(m)).output.weight_MOD{i}];
    end
    datmod(:,2) = num2cell(wemod(:,ffi)*100);
    % Compute Expected Ranking for modalities
    %take 0 columns out of the computation
    wem = wemod(:,1:handles.nfold);
    mw = min(wem);
    mxw = max(wem);
    id0 = find(mw==mxw);
    id0 = id0(mw(id0)==0);
    id = 1:size(wem,2);
    idtr = setdiff(id,id0);
    if isempty(idtr)   % if a modality has always 0 weights
       erwn2=NaN*ones(length(handles.nmods),1);
    else
       wem = wem(:,idtr);
       %deal with NaNs
       [d1,d2]=sort(wem,1,'descend');
       isn=find(isnan(wem(:,1)));
       d3=1:length(isn);
       d4=length(isn)+1:size(d1,1);
       ihn=[d2(d4,:);d2(d3,:)];
       [d1,dwn]=sort(ihn);
       erwn2=zeros(length(handles.nmods),1);
       for i=1:size(wem,1)
           for j=1:size(wem,1)
               tmp=length(find(dwn(i,1:end-1)==j));
               erwn2(i)=erwn2(i)+j*tmp;
           end
       end
       erwn2=erwn2/handles.nfold;
    end
   datmod(:,3) = num2cell(erwn2);
   [vald,idwmod] = sort(wemod(:,ffi),'descend');
   handles.sort_mod= idwmod;
else
    datmod = [];
end
handles.datmod = datmod;  


% Fill table and bar graph if needed
if isfield(handles.PRT.model(mi(m)).output,'weight_ROI') &&... % chosen model has ROI or modality weight values
         ~isempty(handles.PRT.model(mi(m)).output.weight_ROI) &&...
        ~isempty(handles.PRT.model(mi(m)).output.weight_ROI{1}) &&...
        isfield(handles,'wmap') && ~isempty(handles.wmap)
    handles.dattable = dat;
    
    dat = dat(handles.sort_roi,:);
    set(handles.ROItable,'Data',dat);
    set(handles.ROItable,'ColumnEditable',false(1,length(lc)));
    set(handles.ROItable,'ColumnName',lc);
    set(handles.ROItable,'visible','on');
    
    %Bar graph to show decrease in ROI weights
    set(handles.axes1,'visible','on')
    bar(handles.axes1,weights(idwroi));
    set(get(handles.axes1,'XLabel'),'FontWeight','demi')
    set(get(handles.axes1,'XLabel'),'String',xlabel)
    set(get(handles.axes1,'YLabel'),'String',ylabel)
    set(get(handles.axes1,'YLabel'),'FontWeight','demi')
    

   if ~isempty(datmod)
      set(handles.modtable,'Visible','on')
      set(handles.modtable,'Enable','on')
      set(handles.modtable,'Data',datmod(handles.sort_mod,:))
      set(handles.modtable,'ColumnEditable',false(1,size(datmod,2)));
      set(handles.modtable,'ColumnName',{'Modality','Weight (%)','Exp. Ranking'});
   else
      set(handles.modtable,'Visible','off')
      set(handles.modtable,'Enable','off') 
   end    
   
else
    cla(handles.axes1, 'reset');
    set(handles.ROItable,'visible','off');
    set(handles.axes1,'visible','off');
    set(handles.butt_load_labels,'visible','off');
    set(handles.saveweight,'Visible','off')
    set(handles.modtable,'Visible','off')
    set(handles.modtable,'Enable','off')
end
handles.flagmodMKL = flagmodMKL;
handles.model_button = 0;
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function imagemenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imagemenu (see GCBO)
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
child = get(st.handles.weightspanel,'Children');
count = [];
for i = 1:length(child) % Resize colorbar if blobs
    if strcmpi(get(child(i),'Type'),'axes')
        count = [count, i];     
    end
end
if length(count)==4
    pos = get(child(count(1)),'Position');
    set(child(count(1)),'Position',[pos(1)*1.3,pos(2),pos(3),pos(4)*0.9])
end
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
cla(handles.axes1, 'reset');
set(handles.ROItable,'visible','off');
set(handles.axes1,'visible','off');
set(handles.butt_load_labels,'visible','off');
set(handles.saveweight,'visible','off');
set(handles.modtable,'Visible','off')
set(handles.modtable,'Enable','off')
set(handles.loadweight,'String','Load weights map')
set(handles.loadanatomical,'String','Load anatomical img')
handles.noloadw = 0;
gcfchil = get(handles.weightspanel,'Children');
for i=1:length(gcfchil)
    if strcmpi(get(gcfchil(i),'Type'),'axes')
        delete(gcfchil(i)) % Delete any left-over colorbar
    end
end
gcfchil = get(handles.anatomicalpanel,'Children');
for i=1:length(gcfchil)
    if strcmpi(get(gcfchil(i),'Type'),'axes')
        delete(gcfchil(i)) % Delete any left-over colorbar
    end
end
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

print(handles.figure1,ext,[pathname,filesep,b],'-r500')

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

% --- Executes on selection change in table if it is a ROI label.
function disp_weights_CellSelectionCallback(hObject,eventdata,handles)
% hObject    handle to foldmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.selectedcell=eventdata.Indices;
if ~isempty(handles.selectedcell) && handles.selectedcell(2)==1 && ...% ROI label selected
    ~handles.flagmodMKL
    weightbutton_Callback(hObject, eventdata, handles);
else
    handles.selectedcell = [];
end
% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in foldmenu.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to foldmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns foldmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from foldmenu


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in butt_load_labels.
function butt_load_labels_Callback(hObject, eventdata, handles)
% hObject    handle to butt_load_labels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mi = handles.mi;
m = get(handles.classmenu,'Value');
if m==0 % Mac weird error with popup menus
    m=1;
end
mim = get(handles.imagemenu,'Value');
if mim==0
    mim=1;
end
fname=spm_select(1,'.mat','Select the file containing the names of the ROIs');
if isempty(fname)
    disp('No file containing the names of the ROIs found, generic names used')
    handles.labels{mi(m)}{mim}=[];
else
    load(fname)
    try
        handles.labels{mi(m)}{mim}=ROI_names;
    catch
        disp('No variable ROI_names found, generic names used')
    end
end
num_roi = handles.num_roi{mim};
try
    label = handles.labels{mi(m)}{mim}(num_roi);
catch
    if length(num_roi)==length(handles.labels{mi(m)}{mim})  % match between labels and number of regions
        label = handles.labels{mi(m)}{mim}(num_roi);
    elseif length(num_roi)==length(handles.labels{mi(m)}{mim})+1    %missing the 'others' region
        handles.labels{mi(m)}{mim} = [{'others'};handles.labels{mi(m)}{mim}];
        label = handles.labels{mi(m)}{mim}(num_roi);
    else
        warning('prt_ui_disp_weights:LabelsDoNotCorrespondtoROI',... % could not use file
            'Number of labels in .mat does not correspond to number of ROIs');
        label = [];
    end
end
if ~isempty(handles.labels{mi(m)}{mim}) && ~isempty(label)
    dat = handles.dattable;
    dat(:,1) = label;
    handles.dattable = dat;
    dat = dat(handles.sort_roi,:);
    set(handles.ROItable,'Data',dat);
    set(handles.ROItable,'visible','on');
    guidata(hObject, handles);
end





% --- Executes on button press in disp_voxels.
function disp_voxels_Callback(hObject, eventdata, handles)
% hObject    handle to disp_voxels (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of disp_voxels

% Get folds and model indexes
m  = get(handles.classmenu,'Value');
if m==0
    m=1;
end
mi = handles.mi;
handles.nfold = length(handles.PRT.model(mi(m)).output.fold);
ffi = get(handles.foldmenu,'Value')-1;
if ffi == -1 % for Mac issues with popup menus
    ffi = 0;
end
if ffi==0 % for average
    ffi=handles.nfold+1;
end

% Select which image to display: weights per voxels or regions
disp_vox = get(handles.disp_voxels,'Value');
if disp_vox
    fntl = handles.PRT.model(mi(m)).output.weight_img{handles.class};
    set(handles.disp_regions,'Value',0);
else
    fntl = ['ROI_',handles.PRT.model(mi(m)).output.weight_img{handles.class}];%get only first image for now
    set(handles.disp_regions,'Value',1);
    set(handles.disp_voxels,'Value',0);
end

%Load corresponding weight image if it exists
if exist([handles.prtdir,filesep,fntl,'.img'],'file') || ...
        exist([handles.prtdir,filesep,fntl],'file')
    [a,b] = fileparts(fntl);
    handles.wmap = [handles.prtdir,filesep,b,'.img'];
    V               = spm_vol(handles.wmap);
    handles.vols{1} = V;
    handles.noloadw = 1;
    handles.model_button = 0;
    handles.selectedcell = [];
    % Update handles structure
    guidata(hObject, handles);
    weightbutton_Callback(hObject, eventdata, handles);
else
    beep
    disp('Weight image does not exist for this model')
end

    % Update handles structure
    guidata(hObject, handles);

% --- Executes on button press in disp_regions.
function disp_regions_Callback(hObject, eventdata, handles)
% hObject    handle to disp_regions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of disp_regions
% Get folds and model indexes
m  = get(handles.classmenu,'Value');
if m==0
    m=1;
end
mi = handles.mi;
handles.nfold = length(handles.PRT.model(mi(m)).output.fold);
ffi = get(handles.foldmenu,'Value')-1;
if ffi == -1 % for Mac issues with popup menus
    ffi = 0;
end
if ffi==0 % for average
    ffi=handles.nfold+1;
end

% Select which image to display: weights per voxels or regions
disp_reg = get(handles.disp_regions,'Value');
if ~disp_reg
    fntl = handles.PRT.model(mi(m)).output.weight_img{handles.class};
    set(handles.disp_regions,'Value',0);
    set(handles.disp_voxels,'Value',1);
else
    fntl = ['ROI_',handles.PRT.model(mi(m)).output.weight_img{handles.class}];%get only first image for now
    set(handles.disp_regions,'Value',1);
    set(handles.disp_voxels,'Value',0);
end

%Load corresponding weight image if it exists
if exist([handles.prtdir,filesep,fntl,'.img'],'file') || ...
        exist([handles.prtdir,filesep,fntl],'file')
    [a,b] = fileparts(fntl);
    handles.wmap = [handles.prtdir,filesep,b,'.img'];
    V               = spm_vol(handles.wmap);
    handles.vols{1} = V;
    handles.noloadw = 1;
    handles.model_button = 0;
    handles.selectedcell = [];
    % Update handles structure
    guidata(hObject, handles);
    weightbutton_Callback(hObject, eventdata, handles);
else
    beep
    disp('Weight image does not exist for this model')
end

% Update handles structure
guidata(hObject, handles);



% --- Executes when entered data in editable cell(s) in ROItable.
function ROItable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to ROItable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in saveweight.
function saveweight_Callback(hObject, eventdata, handles)
% hObject    handle to saveweight (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% chosen model has ROI values or weights per modality : create the table

%If there are weights per region, then save the weightes for different
%regisons
mim = get(handles.imagemenu,'Value');
if mim==0
    mim=1;
end
num_roi = handles.num_roi{mim};

disp('Saving the sorted list of regions to text file.....>>');
modelname = char(strcat(handles.mnames(get(handles.classmenu,'value')),'_RegionList.txt'));
path = fileparts(handles.pathdir);
weightname=fullfile(path,modelname);
disp(weightname);
fid=fopen(weightname,'w');
fprintf(fid,'%23s %20s %20s %23s\r\n','"ROI Label"','"ROI weight (%)"','"ROI size (vox)"','"Exp. Ranking"');
dat = handles.dattable;
dat = dat(handles.sort_roi,:);
for i = 1:size(num_roi)
    cellvalue= dat(i,:);
    fprintf(fid, '%22s %19f %20d %23f\r\n', cellvalue{:});
end
fclose(fid);
disp('Save Done!');

% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
Tag='WeightResults';
F = findall(allchild(0),'Flat','Tag',Tag);
delete(F);
clear handles
clear F
