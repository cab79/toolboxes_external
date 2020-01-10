function varargout = prt_data_modality(varargin)
% PRT_DATA_MODALITY M-file for prt_data_modality.fig
% 
% PRT_DATA_MODALITY, by itself, creates a new PRT_DATA_MODALITY or raises 
% the existing singleton*.
%
% H = PRT_DATA_MODALITY returns the handle to a new PRT_DATA_MODALITY or 
% the handle to the existing singleton*.
%
% PRT_DATA_MODALITY('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_DATA_MODALITY.M with the given input arguments.
%
% PRT_DATA_MODALITY('Property','Value',...) creates a new PRT_DATA_MODALITY
% or raises the existing singleton*.  Starting from the left, property value
% pairs are applied to the GUI before prt_data_modality_OpeningFcn gets called.  
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_data_modality_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Schrouff
% $Id$

% Edit the above text to modify the response to help prt_data_modality

% Last Modified by GUIDE v2.5 23-Feb-2015 14:52:39

% Begin initialization code - DO NOT EDIT


gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_data_modality_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_data_modality_OutputFcn, ...
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


% --- Executes just before prt_data_modality is made visible.
function prt_data_modality_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_data_modality (see VARARGIN)

% Choose default command line output for prt_data_modality
handles.output = hObject;

%if window already exists, just put it as the current figure
Tag='DDmod';
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

    set(handles.figure1,'Name','PRoNTo :: Specify modality')
    set(handles.design_menu,...
            'String',{'Load SPM.mat','Specify design','No design'},...
            'Value',3);

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
    set(handles.figure1,'DefaultTextFontSize',FS*8,...
        'DefaultUicontrolFontSize',FS*8,...
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
        else
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




    handles.mod=[];
    handles.mod.detrend=1;
    handles.mod.design=0;
    handles.mod.scans=[];
    handles.mod.name={};
    handles.mod.covar=[];
    handles.mod.rt_subj=[];
    handles.subj1=0;

    if ~isempty(varargin) && strcmpi(varargin{1},'UserData')
        %Particular options if you select by 'scans'
        if ~isempty(varargin{2}{2}) && strcmpi(varargin{2}{2}.subj_name,'Scans')
            set(handles.design_menu,'Enable','off')
            set(handles.edit_regt,'Enable','on')
            set(handles.edit_regt,'Visible','on')
            set(handles.edit_covar,'Enable','on')
            set(handles.edit_covar,'Visible','on')
            set(handles.text7,'Visible','on') % Covariates
            set(handles.text6,'Visible','on') %Regression targets
        else
            set(handles.design_menu,'Enable','on')
            set(handles.edit_regt,'Enable','off')
            set(handles.edit_regt,'Visible','off')
            set(handles.edit_covar,'Enable','off')
            set(handles.edit_covar,'Visible','off')
            set(handles.text7,'Visible','off')
            set(handles.text6,'Visible','off')
        end

        if ~isempty(varargin{2}{2}) && isfield(varargin{2}{2},'modality') && ...
                ~isempty(varargin{2}{2}.modality)
            handles.subjmod={varargin{2}{2}.modality(:).mod_name};
            if length(varargin{2})>=3 && ~isempty(varargin{2}{3})
                nlist=varargin{2}{1};
                modsel=varargin{2}{2}.modality(varargin{2}{3});
                valsel=find(strcmpi(modsel.mod_name,nlist));
                set(handles.modname,'String',nlist);
                set(handles.modname,'Value',valsel);
                if isfield(modsel,'detrend')
                    handles.mod.detrend=modsel.detrend;
                else
                    handles.mod.detrend=0;
                end
                handles.mod.design=modsel.design;
                if ~isempty(modsel.design)
                    set(handles.design_menu,'Value',2)
                    handles.desnmenu=2;
                end
                handles.mod.scans=modsel.scans;
                handles.mod.name=modsel.mod_name;
                handles.mod.covar=modsel.covar;
                if ~isempty(modsel.covar)
                    set(handles.edit_covar,'Enable','on')
                    set(handles.edit_covar,'Visible','on')
                    set(handles.text7,'Visible','on')
                    set(handles.edit_covar,'String',num2str(modsel.covar(:)))
                end
                handles.mod.rt_subj=modsel.rt_subj;
                if ~isempty(modsel.rt_subj)
                    set(handles.edit_regt,'Enable','on')
                    set(handles.edit_regt,'Visible','on')
                    set(handles.text6,'Visible','on')
                    set(handles.edit_regt,'String',num2str(modsel.rt_subj(:)))
                end
            else
                nlist=[varargin{2}{1}, {'Enter new'}];
                set(handles.modname,'String',nlist,'Value',length(nlist));  
            end
        else
            nlist=[varargin{2}{1}, {'Enter new'}];
            set(handles.modname,'String',nlist,'Value',length(nlist));  
            handles.subjmod={};
        end
    else
        nlist={'Enter new'};
        handles.subjmod={};
        set(handles.modname,'String',nlist,'Value',1);
    end
    if length(varargin{2})>=4 && ~isempty(varargin{2}{4})
        handles.subj1=varargin{2}{4};
    end
    if length(varargin{2})==5 && ~isempty(varargin{2}{5})
        handles.PRT=varargin{2}{5};
    end
    warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
end
% Update handles structure
guidata(hObject, handles);

% %UIWAIT makes prt_data_modality wait for user response (see UIRESUME)
uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_data_modality_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in modname.
function modname_Callback(hObject, eventdata, handles)
% hObject    handle to modname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns modname contents as cell array
%        contents{get(hObject,'Value')} returns selected item from modname
warning('off','MATLAB:hg:uicontrol:ParameterValuesMustBeValid')
list=get(handles.modname,'String');
%handle particular bug with matlab(7.10.0499) and mac (OSX 10.6.4)
if length(list)==1
    set(handles.modname,'Value',1)
end

if any(strfind(list{get(handles.modname,'Value')}, 'Enter'))
    modname=prt_text_input('Title','Enter modality name');
    if isnumeric(modname)
        return
    end
    if ~any(strcmpi(list,modname))
        nlist=[list;{modname}];
    else
        beep
        disp('This modality has already been set for the selected subject')
        set(handles.modname,'String',list);
        valall=strfind(list,'Enter');
        for i=1:length(valall)
            if ~isempty(valall{i})
                val=i;
                break
            end
        end
        set(handles.modname,'Value',val);
        return               
    end     
    set(handles.modname,'String',nlist);
    set(handles.modname,'Value',length(nlist));
else
    modname=list{get(handles.modname,'Value')};
    if ~isempty(handles.subjmod)
        if any(strcmpi(handles.subjmod,modname))
            set(handles.modname,'String',list);
            valall=strfind(list,'Enter');
            if length(valall)==1 && isempty(valall{1}) %case of reviewing the modality
                val=1;
                set(handles.modname,'Value',val);
            else    %user tries to give a name he already used
                beep
                disp('This modality has already been set for the selected subject')            
                for i=1:length(valall)
                    if ~isempty(valall{i})
                        val=i;
                        break
                    end
                end
                set(handles.modname,'Value',val);
                return
            end
        end
    end         
end

%if a subject was previously entered, then propose to replicate its design
%if the modality is the same
if isstruct(handles.subj1)
    if any(strcmpi(modname, {handles.subj1(:).mod_name}))
        handles.indmods1=find(strcmpi(modname, {handles.subj1(:).mod_name}));
        list=get(handles.design_menu,'String');
        list=[list;{'Replicate design of subject 1'}];
        set(handles.design_menu,'String',list);
    end
end
handles.mod.name=modname;
% Update handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function modname_CreateFcn(hObject, eventdata, handles)
% hObject    handle to modname (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in design_menu.
function design_menu_Callback(hObject, eventdata, handles)
% hObject    handle to design_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns design_menu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from design_menu
choice=get(handles.design_menu,'Value');
%Mac and Matlab versions strange things with popup menus
if choice==0
    choice=handles.desnmenu;
    set(handles.design_menu,'Value')
end
if choice==1
    desn=spm_select(1,'mat','Select SPM.mat file',[],[],'SPM.mat');
    try
        load(desn);
    catch
        beep
        disp('Can not load SPM.mat file');
        return
    end
    conds=struct();
    list={};
    for sspm=1:length(SPM.Sess)
        ncond = length(SPM.Sess(sspm).U);
        for c = 1:ncond
            if sspm==1  %set the conditions
                conds(c).cond_name = SPM.Sess(sspm).U(c).name{1};
                list=[list, {conds(c).cond_name}];
                conds(c).onsets = SPM.Sess(sspm).U(c).ons;
                if length(SPM.Sess(sspm).U(c).dur)==1
                    conds(c).durations=repmat(SPM.Sess(sspm).U(c).dur,...
                        numel(conds(c).onsets),1);
                else
                    conds(c).durations = SPM.Sess(sspm).U(c).dur;
                end
                indcond=ncond;
            else  %for other sessions, look if the conditions are the same
                name_cond=SPM.Sess(sspm).U(c).name{1};
                itoadd=find(strcmpi(name_cond,list));
                if isempty(itoadd)  % if not, then add condition
                    indcond=indcond+1;
                    conds(indcond).cond_name = SPM.Sess(sspm).U(c).name{1};
                    conds(indcond).onsets = SPM.Sess(sspm).U(c).ons;
                    if length(SPM.Sess(sspm).U(c).dur)==1
                        conds(indcond).durations=repmat(SPM.Sess(sspm).U(c).dur,...
                            numel(conds(indcond).onsets),1);
                    else
                        conds(indcond).durations = SPM.Sess(sspm).U(c).dur;
                    end
                    list=[list, {name_cond}];
                else               %if yes, then add the onsets
                    conds(itoadd).onsets    = [conds(itoadd).onsets;SPM.Sess(sspm).U(c).ons];
                    if length(SPM.Sess(sspm).U(c).dur)==1
                        tmp=repmat(SPM.Sess(sspm).U(c).dur,...
                            numel(SPM.Sess(sspm).U(c).ons),1);
                    else
                        tmp = SPM.Sess(sspm).U(c).dur;
                    end
                    conds(itoadd).durations = [conds(itoadd).durations;tmp];
                end
            end
        end
    end    
    if strcmpi(SPM.xBF.UNITS,'scans')
        units=0;
    else
        units=1;
    end
    def=prt_get_defaults('datad');
    if isfield(handles.PRT.group,'hrfoverlap')
        overl=handles.PRT.group.hrfoverlap;
    else
        overl=def.hrfw;
    end
    if isfield(handles.PRT.group,'hrfdelay')
        del=handles.PRT.group.hrfdelay;
    else
        del=def.hrfd;
    end
    %check that TR is the same for all sessions
    if length(unique([SPM.xX.K(:).RT])) ~=1
        beep
        disp('Differents TR found in SPM.mat, please separate sessions')
        return
    end
    desn=prt_check_design(conds,SPM.xX.K(1).RT,units,overl,del);
    desn.covar = [];
elseif choice==2
    if isstruct(handles.mod.design)
        desn=prt_data_conditions('UserData',{handles.mod.design,handles.PRT});
    else
        desn=prt_data_conditions;
    end
elseif choice ==3
    desn=[];
elseif choice==4
    desn=handles.subj1(handles.indmods1).design;
end
handles.mod.design=desn;
if isfield(desn,'covar') && ~isempty(desn.covar)
    set(handles.edit_covar,'String','Entered');
    set(handles.edit_covar,'Visible','on');
    set(handles.text7, 'Visible','on')
end
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function design_menu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to design_menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in getfiles.
function getfiles_Callback(hObject, eventdata, handles)
% hObject    handle to getfiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(handles.mod.scans)
    sel=cellstr(handles.mod.scans);
else
    sel=[];
end
t=spm_select([1 Inf],'image','Select files for the modality',sel);
handles.mod.scans=t;
% Update handles structure
guidata(hObject, handles);


function edit_regt_Callback(hObject, eventdata, handles)
% hObject    handle to edit_regt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_regt as text
%        str2double(get(hObject,'String')) returns contents of edit_regt as a double

%first option of loading: writing the values
rt=get(handles.edit_regt,'String');
if isempty(rt)
    return
end
if strcmp(rt(1),'-'), nrt = 2; else nrt = 1; end
try
    eval(['rte=[',rt,'];'])
catch
    try
        load(char(rt));
    catch
        beep
        disp('Could not load file or read the regression targets')
        disp('Please enter either a .mat file name or enter the values')
        return
    end
    if ~exist('rt_subj','var')
        beep
        sprintf('Regression targets file must contain ''rt_subj'' variable! ')
        disp('Please correct!')
        return
    else
        rte=rt_subj;
    end
end
handles.mod.rt_subj=rte;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_regt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_regt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function edit_covar_Callback(hObject, eventdata, handles)
% hObject    handle to edit_covar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_covar as text
%        str2double(get(hObject,'String')) returns contents of edit_covar as a double
%first option of loading: writing the values
rt=get(handles.edit_covar,'String');
if isempty(rt)
    return
end
try
    eval(['rte=[',rt,'];'])
catch
   try
        load(char(rt));
    catch
        beep
        disp('Could not load file or read the covariate values')
        disp('Please enter either a .mat file name or enter the values')
        return
    end
    if ~exist('R','var')
        beep
        sprintf('Covariates file must contain ''R'' variable! ')
        disp('Please correct!')
        return
    else
        rte=R;
    end
end

handles.mod.covar=rte;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function edit_covar_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_covar (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in okbutton.
function okbutton_Callback(hObject, eventdata, handles)
% hObject    handle to okbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%check that the files were selected
if isempty(handles.mod.scans)
    beep
    disp('Please, select the files before returning')
    return
end

%check that the regression targets have the same number of elements as the
%number of scans
if ~isempty(handles.mod.rt_subj)
    szrt=length(handles.mod.rt_subj);
    if  size(handles.mod.scans,1)~=szrt
        beep
        disp('Number of regression targets must be the number of files selected! ')
        disp('Please correct!')
        return
    end
end
%check that the covariates have the same number of elements as the
%number of scans
if ~isempty(handles.mod.covar)
    szrt=size(handles.mod.covar);
    nsc = size(handles.mod.scans,1);
    ins = find(szrt == nsc);
    if  isempty(ins)
        beep
        disp('Number of covariates must be the number of files selected! ')
        disp('Please correct!')
        return
    else
        if ins~=1 %not the first dimension
            handles.mod.covar = handles.mod.covar';
        end
    end
end
        
modprop=handles.mod;
handles.output=modprop;
% Update handles structure
guidata(hObject, handles);
uiresume(handles.figure1);


% --- Executes on button press in cancelbutton.
function cancelbutton_Callback(hObject, eventdata, handles)
% hObject    handle to cancelbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
uiresume(handles.figure1);


% --- Executes during object creation, after setting all properties.
function text7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
