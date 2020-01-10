function varargout = prt_ui_select_class(varargin)
% PRT_UI_SELECT_CLASS M-file for prt_ui_select_class.fig
% 
% PRT_UI_SELECT_CLASS, by itself, creates a new PRT_UI_SELECT_CLASS or 
% raises the existing singleton*.
%
% H = PRT_UI_SELECT_CLASS returns the handle to a new PRT_UI_SELECT_CLASS 
% or the handle to the existing singleton*.
%
% PRT_UI_SELECT_CLASS('CALLBACK',hObject,eventData,handles,...) calls the 
% local function named CALLBACK in PRT_UI_SELECT_CLASS.M with the given 
% input arguments.
%
% PRT_UI_SELECT_CLASS('Property','Value',...) creates a new PRT_UI_SELECT_CLASS
% or raises the existing singleton*.  Starting from the left, property 
% value pairs are applied to the GUI before prt_ui_select_class_OpeningFcn 
% gets called.  An unrecognized property name or invalid value makes 
% property application stop.  All inputs are passed to prt_ui_select_class_OpeningFcn
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

% Edit the above text to modify the response to help prt_ui_select_class

% Last Modified by GUIDE v2.5 01-Nov-2011 16:35:44

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_select_class_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_select_class_OutputFcn, ...
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


% --- Executes just before prt_ui_select_class is made visible.
function prt_ui_select_class_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_select_class (see VARARGIN)

set(handles.figure1,'Name','PRoNTo :: Specify classes')
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

%set the different fields to disable (will be enabled when choosing the
%number of classes)
set(handles.group_list,'Enable','off')
set(handles.uns_list,'Enable','off')
set(handles.sel_list,'Enable','off')
set(handles.uns_cond_list,'Enable','off')
set(handles.sel_cond_list,'Enable','off')
set(handles.sel_all,'Enable','off')
set(handles.sel_cond_all,'Enable','off');

%get information from the PRT.mat
if ~isempty(varargin) && strcmpi(varargin{1},'UserData')
    handles.dat=varargin{2}{1};
    %get names of groups
    list={handles.dat.group(:).gr_name};
    ng=length(list);
    set(handles.group_list,'String',list)
    
    
    %get the conditions which are common to all groups and subjects for the
    %different modalities comprised in the selected feature set
    indfs=varargin{2}{2};
    handles.indfs=indfs;
    nm=length(handles.dat.fs(indfs).modality);
    modnam={handles.dat.masks(:).mod_name};
    handles.condm=cell(1,3);
    handles.flagcond=1;
    %for each modality, get the conditions which are common to all subjects
    for i=1:nm
        handles.condm{1,1}=get(handles.group_list,'String');
        handles.condm{1,3}=cell(length(get(handles.group_list,'String')),1);
        flag=1;
        for j=1:ng
            handles.condm{1,3}{j}={handles.dat.group(j).subject(:).subj_name};
            for k=1:length(handles.dat.group(j).subject)
                m2= strcmpi(handles.dat.fs(indfs).modality(nm).mod_name,modnam);
                des=handles.dat.group(j).subject(k).modality(m2).design;
                if isstruct(des) && flag
                    if k==1 && j == 1 % [afm] && nm==1
                        lcond={des.conds(:).cond_name};
                    else
                        tocmp={des.conds(:).cond_name};
                        lcond=intersect(lower(lcond),lower(tocmp));
                    end
                else
                    flag=0;
                    lcond={};
                end
            end
        end
    end
    handles.condm{1,2}=lcond;
    if isempty(handles.condm{1,2})
        disp('No conditions were common to all subjects, all groups, all modalities within the feature set')
        disp('Classifying subjects only')
        handles.flagcond=0;
    end
    if length(varargin{2})==2
        handles.flagrev=0;
        set(handles.group_list,'Value',1)
        list={handles.dat.group(1).subject(:).subj_name};
        set(handles.uns_list,'String',list);
        set(handles.sel_list,'String',{});
        % Update handles structure
        guidata(hObject, handles);
        % UIWAIT makes prt_ui_select_class wait for user response (see UIRESUME)
        uiwait(handles.figure1);
    else
        %review the model
        handles.flagrev=1;
        indm=varargin{2}{3};
        mod=handles.dat.model(indm).input;
        %classification case
        nc=length(mod.class);
        set(handles.num_class,'String',num2str(nc))
        set(handles.num_class,'Enable','off')
        set(handles.edit2,'Enable','off')
        handles.clas=cell(nc,4);
        %set the names of the classes in the pop_class
        cl={};
        for i=1:nc
            cl=[cl,{['Class ',num2str(i)]}];
        end
        set(handles.pop_class,'String',cl);
        set(handles.pop_class,'Value',1);
        for i=1:nc
            handles.class(i).class_name=mod.class(i).class_name;
            handles.clas{i,1}=1:length(handles.condm{1,1});
            handles.clas{i,2}=cell(length(handles.condm{1,1}),2);
            list=get(handles.group_list,'String');
            for j=1:length(list)
                indg=find(strcmpi(list{j},{mod.class(i).group(:).gr_name}));
                if~isempty(indg)
                    sel=[mod.class(i).group(indg).subj(:).num];
                    all=1:length(handles.condm{1,3}{j});
                    handles.clas{i,2}{j,1}=setdiff(all,sel);
                    handles.clas{i,2}{j,2}=sel;
                else
                    handles.clas{i,2}{j,1}=1:length(handles.condm{1,3}{j});
                    handles.clas{i,2}{j,2}=0;
                end
                if isempty(handles.clas{i,2}{j,1})
                    handles.clas{i,2}{j,1}=0;
                end
                if isempty(handles.clas{i,2}{j,2})
                    handles.clas{i,2}{j,2}=0;
                end
            end
            if handles.flagcond
                listc=handles.condm{1,2};
                selc=[];
                if isfield(mod.class(i).group(1).subj(1).modality(1),'all_cond')
                    selc=1:length(listc);
                elseif isfield(mod.class(i).group(1).subj(1).modality(1),'conds')
                    for icc=1:length(listc)
                        indcc=strcmpi(listc{icc},{mod.class(i).group(1).subj(1).modality(1).conds(:).cond_name});
                        if any(indcc)
                            selc=[selc,icc];
                        end
                    end
                end
                allc=1:length(handles.condm{1,2});
                handles.clas{i,3}=setdiff(allc,selc);
                handles.clas{i,4}=selc;
                if isempty(handles.clas{i,3})
                    handles.clas{i,3}=0;
                end
                if isempty(handles.clas{i,4})
                    handles.clas{i,3}=0;
                end
            end
            handles.class(i).class_name=mod.class(i).class_name;            
        end
        set(handles.group_list,'Enable','on');
        % Update handles structure
        guidata(hObject, handles);
    end
end






% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_select_class_OutputFcn(hObject, eventdata, handles) 
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
if isfield(handles,'figure1') && ~handles.flagrev
    delete(handles.figure1)
end



function num_class_Callback(hObject, eventdata, handles)
% hObject    handle to num_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of num_class as text
%        str2double(get(hObject,'String')) returns contents of num_class as a double
ncl=str2double(get(handles.num_class,'String'));
if isnan(ncl)
    return
end
%set the names of the classes in the pop_class
cl={};
for i=1:ncl
    cl=[cl,{['Class ',num2str(i)]}];
end
set(handles.pop_class,'String',cl);
set(handles.pop_class,'Value',1);
%For each class, build a cell containing the indexes of
%the selected groups and conditions
handles.clas=cell(ncl,4);
for i=1:ncl
    handles.clas{i,1}=1:length(handles.condm{1,1});
    handles.clas{i,2}=cell(length(handles.condm{1,1}),2);
    for j=1:length(get(handles.group_list,'String'))
        handles.clas{i,2}{j,1}=1:length(handles.condm{1,3}{j});
        handles.clas{i,2}{j,2}=0;
    end
    if handles.flagcond
        handles.clas{i,3}=1:length(handles.condm{1,2});
        handles.clas{i,4}=0;
    end
    handles.class(i).class_name=cl{i};
end
set(handles.edit2,'String',handles.class(1).class_name)
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function num_class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%change the name of the classes
function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double
indc=get(handles.pop_class,'Value');
handles.class(indc).class_name=get(handles.edit2,'String');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in pop_class.
function pop_class_Callback(hObject, eventdata, handles)
% hObject    handle to pop_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pop_class contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_class
vc=get(handles.pop_class,'Value');
if vc==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_class,'Value',1)
    vc=1;
end
cg=get(handles.group_list,'Value');
list=handles.condm{1,3}{cg};
clist=handles.condm{1,2};
set(handles.edit2,'String',handles.class(vc).class_name)
%set subjects lists
if handles.clas{vc,2}{cg,1}~=0
    set(handles.uns_list,'Value',1);
    set(handles.uns_list,'String',list(handles.clas{vc,2}{cg,1}));
else
    set(handles.uns_list,'Value',0);
    set(handles.uns_list,'String',{});
end
if handles.clas{vc,2}{cg,2}~=0
    set(handles.sel_list,'Value',1);
    set(handles.sel_list,'String',list(handles.clas{vc,2}{cg,2}));
else
    set(handles.sel_list,'Value',0);
    set(handles.sel_list,'String',{});
end
%set conditions lists
if handles.flagcond
    if handles.clas{vc,3}~=0
        set(handles.uns_cond_list,'Value',1);
        set(handles.uns_cond_list,'String',clist(handles.clas{vc,3}));
    else
        set(handles.uns_cond_list,'Value',0);
        set(handles.uns_cond_list,'String',{});
    end
    if handles.clas{vc,4}~=0
        set(handles.sel_cond_list,'Value',1);
        set(handles.sel_cond_list,'String',clist(handles.clas{vc,4}));
    else
        set(handles.sel_cond_list,'Value',0);
        set(handles.sel_cond_list,'String',{});
    end
    if ~handles.flagrev
        set(handles.uns_cond_list,'Enable','on');
        set(handles.sel_cond_list,'Enable','on');
        set(handles.sel_cond_all,'Enable','on');
    end
end
if ~handles.flagrev
    set(handles.uns_list,'Enable','on');
    set(handles.sel_list,'Enable','on');
    set(handles.sel_all,'Enable','on');
end
set(handles.group_list,'Enable','on');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function pop_class_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_class (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in group_list.
function group_list_Callback(hObject, eventdata, handles)
% hObject    handle to group_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns group_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from group_list
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
cg=get(handles.group_list,'Value');
list=handles.condm{1,3}{cg};
%set subjects lists
if handles.clas{cl,2}{cg,1}~=0
    set(handles.uns_list,'Value',1);
    set(handles.uns_list,'String',list(handles.clas{cl,2}{cg,1}));
else
    set(handles.uns_list,'Value',0);
    set(handles.uns_list,'String',{});
end
if handles.clas{cl,2}{cg,2}~=0
    set(handles.sel_list,'Value',1);
    set(handles.sel_list,'String',list(handles.clas{cl,2}{cg,2}));
else
    set(handles.sel_list,'Value',0);
    set(handles.sel_list,'String',{});
end
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

% --- Executes on selection change in uns_list.
function uns_list_Callback(hObject, eventdata, handles)
% hObject    handle to uns_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns uns_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from uns_list
val=get(handles.uns_list,'Value');
induns=1:length(get(handles.uns_list,'String'));
indok=setdiff(induns,val);
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
cg=get(handles.group_list,'Value');
if handles.clas{cl,2}{cg,2}==0
    handles.clas{cl,2}{cg,2}=handles.clas{cl,2}{cg,1}(val);
else
    handles.clas{cl,2}{cg,2}=[handles.clas{cl,2}{cg,2}, handles.clas{cl,2}{cg,1}(val)];
end
if isempty(indok)
    handles.clas{cl,2}{cg,1}=0;
else
    handles.clas{cl,2}{cg,1}=handles.clas{cl,2}{cg,1}(indok);
end
list=handles.condm{1,3}{cg};
%set subjects lists
if handles.clas{cl,2}{cg,1}~=0
    set(handles.uns_list,'Value',1);
    set(handles.uns_list,'String',list(handles.clas{cl,2}{cg,1}));
else
    set(handles.uns_list,'Value',0);
    set(handles.uns_list,'String',{});
end
if handles.clas{cl,2}{cg,2}~=0
    set(handles.sel_list,'Value',1);
    set(handles.sel_list,'String',list(handles.clas{cl,2}{cg,2}));
else
    set(handles.sel_list,'Value',0);
    set(handles.sel_list,'String',{});
end

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uns_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uns_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sel_list.
function sel_list_Callback(hObject, eventdata, handles)
% hObject    handle to sel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sel_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_list
val=get(handles.sel_list,'Value');
indsel=1:length(get(handles.sel_list,'String'));
indok=setdiff(indsel,val);
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
cg=get(handles.group_list,'Value');
if handles.clas{cl,2}{cg,1}==0
    handles.clas{cl,2}{cg,1}=handles.clas{cl,2}{cg,2}(val);
else
    handles.clas{cl,2}{cg,1}=[handles.clas{cl,2}{cg,1}, handles.clas{cl,2}{cg,2}(val)];
end
if isempty(indok)
    handles.clas{cl,2}{cg,2}=0;
else
    handles.clas{cl,2}{cg,2}=handles.clas{cl,2}{cg,2}(indok);
end
list=handles.condm{1,3}{cg};
%set subjects lists
if handles.clas{cl,2}{cg,1}~=0
    set(handles.uns_list,'String',list(handles.clas{cl,2}{cg,1}));
    set(handles.uns_list,'Value',length(get(handles.uns_list,'String')));
else
    set(handles.uns_list,'Value',0);
    set(handles.uns_list,'String',{});
end
if handles.clas{cl,2}{cg,2}~=0
    set(handles.sel_list,'String',list(handles.clas{cl,2}{cg,2}));
    set(handles.sel_list,'Value',length(get(handles.sel_list,'String')));
else
    set(handles.sel_list,'Value',0);
    set(handles.sel_list,'String',{});
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sel_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sel_all.
function sel_all_Callback(hObject, eventdata, handles)
% hObject    handle to sel_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
cg=get(handles.group_list,'Value');
list=handles.condm{1,3}{cg,1};
indsel=1:length(list);
handles.clas{cl,2}{cg,2}=indsel;
handles.clas{cl,2}{cg,1}=0;
set(handles.uns_list,'String',{});
set(handles.uns_list,'Value',0);
set(handles.sel_list,'String',list(handles.clas{cl,2}{cg,2}));
set(handles.sel_list,'Value',length(get(handles.sel_list,'String')));
% Update handles structure
guidata(hObject, handles);


% --- Executes on selection change in uns_cond_list.
function uns_cond_list_Callback(hObject, eventdata, handles)
% hObject    handle to uns_cond_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns uns_cond_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from uns_cond_list
val=get(handles.uns_cond_list,'Value');
induns=1:length(get(handles.uns_cond_list,'String'));
indok=setdiff(induns,val);
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
if handles.clas{cl,4}==0
    handles.clas{cl,4}=handles.clas{cl,3}(val);
else
    handles.clas{cl,4}=[handles.clas{cl,4}, handles.clas{cl,3}(val)];
end
if isempty(indok)
    handles.clas{cl,3}=0;
else
    handles.clas{cl,3}=handles.clas{cl,3}(indok);
end
%set conditions lists
clist=handles.condm{1,2};
if handles.clas{cl,3}~=0
    set(handles.uns_cond_list,'Value',1);
    set(handles.uns_cond_list,'String',clist(handles.clas{cl,3}));
else
    set(handles.uns_cond_list,'Value',0);
    set(handles.uns_cond_list,'String',{});
end
if handles.clas{cl,4}~=0
    set(handles.sel_cond_list,'Value',1);
    set(handles.sel_cond_list,'String',clist(handles.clas{cl,4}));
else
    set(handles.sel_cond_list,'Value',0);
    set(handles.sel_cond_list,'String',{});
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function uns_cond_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uns_cond_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in sel_cond_list.
function sel_cond_list_Callback(hObject, eventdata, handles)
% hObject    handle to sel_cond_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns sel_cond_list contents as cell array
%        contents{get(hObject,'Value')} returns selected item from sel_cond_list
val=get(handles.sel_cond_list,'Value');
induns=1:length(get(handles.sel_cond_list,'String'));
indok=setdiff(induns,val);
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
if handles.clas{cl,3}==0
    handles.clas{cl,3}=handles.clas{cl,4}(val);
else
    handles.clas{cl,3}=[handles.clas{cl,3}, handles.clas{cl,4}(val)];
end
if isempty(indok)
    handles.clas{cl,4}=0;
else
    handles.clas{cl,4}=handles.clas{cl,4}(indok);
end
%set conditions lists
clist=handles.condm{1,2};
if handles.clas{cl,3}~=0
    set(handles.uns_cond_list,'String',clist(handles.clas{cl,3}));
    set(handles.uns_cond_list,'Value',length(get(handles.uns_cond_list,'String')));
else
    set(handles.uns_cond_list,'Value',0);
    set(handles.uns_cond_list,'String',{});
end
if handles.clas{cl,4}~=0
    set(handles.sel_cond_list,'String',clist(handles.clas{cl,4}));
    set(handles.sel_cond_list,'Value',length(get(handles.sel_cond_list,'String')));
else
    set(handles.sel_cond_list,'Value',0);
    set(handles.sel_cond_list,'String',{});
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function sel_cond_list_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sel_cond_list (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sel_cond_all.
function sel_cond_all_Callback(hObject, eventdata, handles)
% hObject    handle to sel_cond_all (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cl=get(handles.pop_class,'Value');
if cl==0
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
    set(handles.pop_classes,'Value',1)
end
cl=get(handles.pop_class,'Value');
list=handles.condm{1,2};
indsel=1:length(list);
handles.clas{cl,4}=indsel;
handles.clas{cl,3}=0;
set(handles.uns_cond_list,'String',{});
set(handles.uns_cond_list,'Value',0);
set(handles.sel_cond_list,'String',list(handles.clas{cl,4}));
set(handles.sel_cond_list,'Value',length(get(handles.sel_cond_list,'String')));
% Update handles structure
guidata(hObject, handles);



% --- Executes on button press in done_button.
function done_button_Callback(hObject, eventdata, handles)
% hObject    handle to done_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.flagrev
    delete(handles.figure1)
    return
end
ncc=[];
scc=zeros(size(handles.clas,1),1);
for i=1:size(handles.clas,1)
    list=get(handles.group_list,'String');
    flag=0;
    ncs=[];
    for g=1:length(list)
        scount=1;
        g2=find(strcmpi(list{g},{handles.dat.group(:).gr_name}));
        sids=handles.clas{i,2}{g,2};
        if ~isempty(sids) && any(sids)
            handles.class(i).group(g2).gr_name=list{g};
            flag=1;
            for s=1:length(sids)
                handles.class(i).group(g2).subj(scount).num=sids(s);
                listm={handles.dat.fs(handles.indfs).modality(:).mod_name};
                for m=1:length(listm)
                    handles.class(i).group(g2).subj(scount).modality(m).mod_name=listm{m};
                    if isempty(handles.condm{1,2})
                        handles.class(i).group(g2).subj(scount).modality(m).all_scans=true;
                        ncs=[ncs;0];
                    else %design and conditions selected
                        if isempty(handles.clas{i,3}) || any(handles.clas{i,3}==0) %all conditions were selected
                            handles.class(i).group(g2).subj(scount).modality(m).all_cond=true;
                        else
                            for ic=1:length(handles.clas{i,4})
                                handles.class(i).group(g2).subj(scount).modality(m).conds(ic).cond_name= ...
                                    handles.condm{1,2}{handles.clas{i,4}(ic)};
                            end
                        end
                        ncs=[ncs;length(handles.clas{i,4})];
                    end
                end
                scount=scount+1;
            end
            scc(i)=scc(i)+scount;
        end
    end
    if ~flag
        beep
        sprintf('No subjects found in the definition of class %d',i)
        disp('Please select subjects (and conditions) for that class')
        return
    end
    if length(unique(ncs))~=1
        beep
        sprintf('Different numbers of conditions found in the definition of class %d',i)
        disp('Please select the same conditions for each subject/group of that class')
        return
    else
        ncc=[ncc;unique(ncs)];
    end
    if ~isempty(find(ncc==0)) && length(find(ncc==0))~=length(ncc)
        beep
        sprintf('Class %d does not have the same number of conditions as class 1',i)
        disp('Please select either at least one condition for each class or none')
        return
    end
end

if ~any(ncc) %no conditions specified
    handles.design=0;
else
    handles.design=1;
end
if length(unique(scc))~=1  %different numbers of subject per class
    handles.loospg=0;  %no leave-one-subject per group-out CV
else
    handles.loospg=1;
end

aa=struct();
aa.class=handles.class;
aa.design=handles.design;
aa.loospg=handles.loospg;
%get names of selected groups and conditions for custom CV GUI
lgi=[];
lci=[];
for i=1:size(handles.clas,1)
    %get which groups
    d=handles.clas{i,2};
    for j=1:size(d,1)
        if ~isempty(d{j,2}) && any(d{j,2}~=0)
            lgi=[lgi,j]; %d{j,2}
        end
    end
    %get which conditions
    d=handles.clas{i,4};
    if ~isempty(d) && any(d~=0)
        lci=[lci,d];
    end
end
lg=handles.condm{1}(lgi);
lc=handles.condm{2}(lci);
legends=struct();
legends.lg=lg;
legends.lgi=lgi;
legends.lci=lci;
legends.lc=lc;
aa.legends=legends;
handles.output=aa;
% Update handles structure
guidata(hObject, handles);
if ~handles.flagrev
    uiresume(handles.figure1)
end
