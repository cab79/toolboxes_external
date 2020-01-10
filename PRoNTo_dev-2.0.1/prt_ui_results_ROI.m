function varargout = prt_ui_results_ROI(varargin)
% PRT_UI_RESULTS_ROI M-file for prt_ui_results_ROI.fig
% 
% PRT_UI_RESULTS_ROI, by itself, creates a new PRT_UI_RESULTS_ROI or 
% raises the existing singleton*.
%
% H = PRT_UI_RESULTS_ROI returns the handle to a new PRT_UI_RESULTS_ROI 
% or the handle to the existing singleton*.
%
% PRT_UI_RESULTS_ROI('CALLBACK',hObject,eventData,handles,...) calls the 
% local function named CALLBACK in PRT_UI_RESULTS_ROI.M with the given 
% input arguments.
%
% PRT_UI_RESULTS_ROI('Property','Value',...) creates a new 
% PRT_UI_RESULTS_ROI or raises the existing singleton*.  Starting from the
% left, property value pairs are applied to the GUI before 
% prt_ui_results_ROI_OpeningFcn gets called.  An unrecognized property name
% or invalid value makes property application stop.  All inputs are passed 
% to prt_ui_results_ROI_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%  instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%_______________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J.Schrouff
% $Id: prt_ui_results_ROI.m 555 2012-06-12 08:44:51Z cphillip $

% Edit the above text to modify the response to help prt_ui_results_ROI

% Last Modified by GUIDE v2.5 18-Feb-2013 10:31:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_results_ROI_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_results_ROI_OutputFcn, ...
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


% --- Executes just before prt_ui_results_ROI is made visible.
function prt_ui_results_ROI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_results_ROI (see VARARGIN)

% Choose default command line output for prt_ui_results_ROI
handles.output = hObject;

%if window already exists, just put it as the current figure
Tag='ResROI';
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
    %build figure when it doesn't exist yet
    set(handles.figure1,'Name','PRoNTo :: ROI weights')
    %set size of the window, taking screen resolution and platform into account
    S0= spm('WinSize','0',1);   %-Screen size (of the current monitor)
    if ispc
        PF ='MS Sans Serif';
    else
        PF = spm_platform('fonts');     %-Font names (for this platform)
        PF = PF.helvetica;
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
                    'radiobutton','checkbox','listbox'})))
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
    
    %fill table with UserData
    if ~isempty(varargin) && strcmpi(varargin{1},'UserData')
        LR=varargin{2}{1};
        pHN=varargin{2}{2};
        drwn=varargin{2}{3};
        erwn=varargin{2}{4};
        SN=varargin{2}{5};
        P_oth=varargin{2}{6};
        oth_w=varargin{2}{7};
    else
        f=spm_select(1,'.mat','Select .mat created from atlas weights');
        if ~isempty(f)
            load(f);
            try
                pHN=NW_roi;
            catch
                error('PRoNTo:prt_ui_results_ROI:Nodata',...
                    'No ROI names or weight values found')
            end
        else
            error('PRoNTo:prt_ui_results_ROI:Nomat',...
                 'No file selected')
        end
        %future: load the saved .mat file
    end
    szn=length(LR);
    dat=cell(szn,4);
    for i=1:szn
        dat{i,1}=LR{i};
        dat{i,2}=pHN(i,end);
        dat{i,3}=erwn(i);
        dat{i,4}=SN(i,end);
    end
    handles.LR=LR;
    handles.pHN=pHN;
    handles.SN=SN;
    handles.rwn=erwn;
    handles.dat=dat;
    set(handles.ROItable,'Data',dat);
    set(handles.ROItable,'ColumnName',{'ROI Name','NW_roi (in %)',...
        'Ranking','Pos weights (in %)'});
    set(handles.ROItable,'ColumnEditable',[false,false,false,false]);
    set(handles.ROItable,'ColumnWidth',{130,80,80,80});
    set(handles.ROItable,'ColumnFormat',{'char','numeric',...
        'numeric','numeric'});
    
    %Set the values for the 'others' region
    handles.P_oth=P_oth;
    handles.oth_w=oth_w;
    set(handles.vol_oth,'String',num2str(P_oth*100))
    set(handles.w_oth,'String',num2str(oth_w))
    
    %get the number of folds to create the popup menu
    handles.nfolds=size(pHN,2)-1;
    list=get(handles.foldp,'String');
    for ii=1:handles.nfolds
        list=[list,{['Fold ',num2str(ii)]}];
    end
    set(handles.foldp,'String',list)
    set(handles.foldp,'Value',1)
    
    % Update Histogram of weights
    set(handles.figure1, 'Currentaxes',handles.axes1)
    H1 = bar(1:size(dat,1),pHN(:,end));
    set(H1,'BarWidth',0.8)
    set(gca,'xlim',[0 size(dat,1)+1])
    title('Histogram of normalized weights per region','FontWeight','bold')
    set(gca,'XTickLabel',{''})
    set(gca,'XTick',0)
    legend('NW_{roi}','Location','NorthEast')
    
    %Update histogram of ranking distances between each fold and the average
    set(handles.figure1, 'Currentaxes',handles.axes2)
    h=bar(handles.axes2,drwn');
    set(get(handles.axes2,'Title'),'String','Ranking distance to average across folds',...
        'Fontweight','bold')
    xlabel(handles.axes2,'Folds','Fontweight','bold')
    set(handles.axes2,'Xtick',1:handles.nfolds)
    set(handles.axes2,'XTickLabel',{''})
end
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_results_ROI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

varargout{1} = handles.output;



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

% Update handles structure
guidata(hObject, handles);

% --- Executes on selection change in foldp.
function foldp_Callback(hObject, eventdata, handles)
% hObject    handle to foldp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns foldp contents as cell array
%        contents{get(hObject,'Value')} returns selected item from foldp
v=get(handles.foldp,'Value');
if v==0 || isempty(v)
    v=1;
end

szn=length(handles.LR);
dat=cell(szn,4);
for i=1:szn
    dat{i,1}=handles.LR{i};
    if v==1 %dispay the average
        dat{i,2}=handles.pHN(i,end);
        dat{i,4}=handles.SN(i,end);
    else
        dat{i,2}=handles.pHN(i,v-1);
        dat{i,4}=handles.SN(i,v-1);
    end
    dat{i,3}=handles.rwn(i);
end
set(handles.ROItable,'Data',dat)

% Update Histogram of weights
set(handles.figure1, 'Currentaxes',handles.axes1)
H1 = bar(1:size(dat,1),[dat{:,2}]);
set(H1,'BarWidth',0.8)
set(gca,'xlim',[0 size(dat,1)+1])
title('Histogram of normalized weights per region','FontWeight','bold')
set(gca,'XTickLabel',{''})
set(gca,'XTick',0)
legend('NW_{roi}','Location','NorthEast')

% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function foldp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to foldp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in sortHN.
function sortHN_Callback(hObject, eventdata, handles)
% hObject    handle to sortHN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
v=get(handles.foldp,'Value');
if v==0 || isempty(v)
    v=1;
end
if v==1
    [d1,d2]=sort(handles.pHN(:,end),'descend');
else
    [d1,d2]=sort(handles.pHN(:,v-1),'descend');
end
isn=find(isnan(handles.pHN(:,1)));
d3=1:length(isn);
d4=length(isn)+1:size(d1,1);
ihn=[d2(d4,:);d2(d3,:)];
szn=length(handles.LR);
dat=cell(szn,4);
for i=1:szn
    dat{i,1}=handles.LR{i};
    if v==1 %dispay the average
        dat{i,2}=handles.pHN(i,end);
        dat{i,4}=handles.SN(i,end);
    else
        dat{i,2}=handles.pHN(i,v-1);
        dat{i,4}=handles.SN(i,v-1);
    end
    dat{i,3}=handles.rwn(i);
end
dat=dat(ihn,:);
set(handles.ROItable,'Data',dat)

% Update Histogram of weights
set(handles.figure1, 'Currentaxes',handles.axes1)
H1 = bar(1:size(dat,1),[dat{:,2}]);
set(H1,'BarWidth',0.8)
set(gca,'xlim',[0 size(dat,1)+1])
title('Histogram of normalized weights per region','FontWeight','bold')
set(gca,'XTickLabel',{''})
set(gca,'XTick',0)
legend('NW_{roi}','Location','NorthEast')

% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in quitbutt.
function quitbutt_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% The figure can be deleted now
if isfield(handles,'figure1')
    delete(handles.figure1);
end



function vol_oth_Callback(hObject, eventdata, handles)
% hObject    handle to vol_oth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of vol_oth as text
%        str2double(get(hObject,'String')) returns contents of vol_oth as a
%        double


% --- Executes during object creation, after setting all properties.
function vol_oth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to vol_oth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function w_oth_Callback(hObject, eventdata, handles)
% hObject    handle to w_oth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of w_oth as text
%        str2double(get(hObject,'String')) returns contents of w_oth as a double


% --- Executes during object creation, after setting all properties.
function w_oth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to w_oth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
