function varargout = prt_ui_reviewCV(varargin)
% PRT_UI_REVIEWCV M-file for prt_ui_reviewCV.fig
% 
% PRT_UI_REVIEWCV, by itself, creates a new PRT_UI_REVIEWCV or raises the 
% existing singleton*.
%
% H = PRT_UI_REVIEWCV returns the handle to a new PRT_UI_REVIEWCV or the 
% handle to the existing singleton*.
%
% PRT_UI_REVIEWCV('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_REVIEWCV.M with the given input 
% arguments.
%
% PRT_UI_REVIEWCV('Property','Value',...) creates a new PRT_UI_REVIEWCV or 
% raises the existing singleton*.  Starting from the left, property value 
% pairs are applied to the GUI before prt_ui_reviewCV_OpeningFcn gets 
% called.  An unrecognized property name or invalid value makes property 
% application stop.  All inputs are passed to prt_ui_reviewCV_OpeningFcn 
% via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by J. Schrouff
% $Id$

% Edit the above text to modify the response to help prt_ui_reviewCV

% Last Modified by GUIDE v2.5 28-Mar-2012 15:59:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prt_ui_reviewCV_OpeningFcn, ...
    'gui_OutputFcn',  @prt_ui_reviewCV_OutputFcn, ...
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


% --- Executes just before prt_ui_reviewCV is made visible.
function prt_ui_reviewCV_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_reviewCV (see VARARGIN)

% Choose default command line output for prt_ui_reviewCV
handles.output = hObject;
%if window already exists, just put it as the current figure
Tag='CVrev';
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
    
    set(handles.figure1,'Name','PRoNTo :: Review Cross-Validation')
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
    % set(handles.figure1,'DefaultTextFontSize',FS*12,...
    %     'DefaultUicontrolFontSize',FS*12,...
    %     'DefaultTextFontName',PF,...
    %     'DefaultAxesFontName',PF,...
    %     'DefaultUicontrolFontName',PF)
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
        if ~strcmpi(get(aa(i),'type'),'uimenu')
            set(aa(i),'FontUnits','pixel')
            xf=get(aa(i),'FontSize');
            set(aa(i),'FontSize',ceil(FS*xf),'FontName',PF,...
                'Units','normalized')
        end
    end


if ~isempty(varargin{1}) && strcmpi(varargin{1},'UserData')
    handles.PRT=varargin{2}{1};
    handles.indm=varargin{2}{2};
    handles.indf=varargin{2}{3};
    handles.prtdir=varargin{2}{4};
else
    beep
    disp('The PRT, index of the model and index of feature set should be entered')
    return
end
end
tcv = handles.PRT.model(handles.indm).input.cv_type;
kcv = handles.PRT.model(handles.indm).input.cv_k;
if kcv>0 %k-folds CV
    cv='k-folds on ';
else
    cv='Leave One ';
end
if strcmpi(tcv,'lobo')
    cv = [cv, 'Block Out'];
elseif strcmpi(tcv,'loro')
    cv = [cv, 'Run Out'];
elseif strcmpi(tcv,'loso')
    cv = [cv, 'Subject Out'];
elseif strcmpi(tcv,'losgo')
    cv = [cv, 'Subject per Group Out'];
elseif strcmpi(tcv,'custom')
    cv = 'Custom';
end
set(handles.selCV,'String',cv);
% Update handles structure
guidata(hObject, handles);
disp_cv(hObject,handles,handles.indm,handles.indf);



% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_reviewCV_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


%--------------------------------------------------------------------------
%-------------------------  Subfunctions ----------------------------------
%--------------------------------------------------------------------------

function disp_cv(hObject,handles,indm,indf)

cla(handles.axes1)
cla(handles.axes2)
cla(handles.axes3)

%Plot the id_mat in the left part of the window
set(handles.figure1,'CurrentAxes',handles.axes1)
dat=handles.PRT.fs(indf).id_mat(:,1:6);
for i=1:6
    dat(:,i)=dat(:,i)./max(dat(:,i));
end
imagesc(dat);
set(gca,'XTick',1:6)
set(gca,'XTickLabel',{'Group','Subject','Modality','Condition','Block','Scans'},...
    'FontWeight','demi','FontSize',9);
set(gca,'YTickLabel',{})
colorbar('Location','WestOutside')
set(get(gca,'Title'),'String','Feature set','FontWeight','bold')

%Plot the CV matrix in the right part of the window
set(handles.figure1,'CurrentAxes',handles.axes2)
CV_mat_full=zeros(size(handles.PRT.fs(indf).id_mat,1),...
    size(handles.PRT.model(indm).input.cv_mat,2));
xticksl=cell(1,size(handles.PRT.model(indm).input.cv_mat,2));
for i=1:size(handles.PRT.model(indm).input.cv_mat,2)
    CV_mat_full(handles.PRT.model(indm).input.samp_idx,i)=handles.PRT.model(indm).input.cv_mat(:,i);
    xticksl{i}=num2str(i);
end
%inverting unused and test for colour purposes
if max(max(CV_mat_full)-min(CV_mat_full))>1
    indun=find(CV_mat_full==0);
    indt=find(CV_mat_full==2);
    CV_mat_full(indun)=2;
    CV_mat_full(indt)=0;
else
    indtr=find(CV_mat_full==1);
    indt=find(CV_mat_full==2);
    CV_mat_full(indtr)=2;
    CV_mat_full(indt)=1;
end

set(gca,'FontWeight','bold')
xlabel('CV Folds','fontweight','demi')
imagesc(CV_mat_full);
set(gca,'XTick',1:i)
set(gca,'YTickLabel',{})
set(gca,'XTickLabel',xticksl,'FontWeight','demi','FontSize',9)
colormap(gray)
set(get(gca,'Title'),'String','Cross-Validation','FontWeight','bold')

%Plot the 'legend' corresponding to the CV matrix in the right bottom part
set(handles.figure1,'CurrentAxes',handles.axes3)
if max(max(CV_mat_full)-min(CV_mat_full))>1
    leg=[0; 1; 2];
    imagesc(leg);
    set(gca,'YTick',[1,2,3])
    set(gca,'YTickLabel',{'Test','Train','Unused'});
else
    leg=[1; 2];
    imagesc(leg);
    set(gca,'YTick',[1,2])
    set(gca,'YTickLabel',{'Test','Train'});
end
set(gca,'YAxisLocation','right')
set(gca,'XTickLabel',{})

% Update handles structure
guidata(hObject, handles);


% --------------------------------------------------------------------
function savemenu_Callback(hObject, eventdata, handles)
% hObject    handle to savemenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
wd=cd;
cd(handles.prtdir)
modname=handles.PRT.model(handles.indm).model_name;
[filename, pathname] = uiputfile( ...
{'*.png','Portable Network Graphics (*.png)';...
 '*.jpeg','JPEG figure (*.jpeg)';...
 '*.tiff','Compressed TIFF figure (*.tiff)';... 
 '*.fig','Matlab figure (*.fig)';...
 '*.pdf','Color PDF file (*.pdf)';...
 '*.epsc',  'Encapsulated PostScript (*.eps)'},...
 'Save figure as',['Cross_Validation_',modname,'.png']);
[a,b,c]=fileparts(filename);
ext=['-d',c(2:end)];

% Set the color of the different backgrounds and figure parameters to white
cf=get(handles.figure1,'Color');
set(handles.figure1,'Color',[1,1,1])
aa=get(handles.figure1,'children');
c=zeros(length(aa),3);
cb=cell(length(aa));
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        bb=get(aa(i),'children');
        cb{i}=zeros(length(bb),3);
        if ~isempty(bb)
            for j=1:length(bb)
                try
                    cb{i}(j,:)=get(bb(j),'BackgroundColor');
                    set(bb(j),'BackgroundColor',[1 1 1]);
                end
            end
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            c(i,:)=get(aa(i),'BackgroundColor');
            set(aa(i),'BackgroundColor',[1 1 1]);
        end
    end
end

print(handles.figure1,ext,[pathname,filesep,b],'-r500')

% Set the color of the different backgrounds and figure parameters to white
set(handles.figure1,'Color',cf)
aa=get(handles.figure1,'children');
for i=1:length(aa)
    if strcmpi(get(aa(i),'type'),'uipanel')
        bb=get(aa(i),'children');
        if ~isempty(bb)
            for j=1:length(bb)
                set(bb(j),'BackgroundColor',cb{i}(j,:));
            end
        end
    end
    if ~strcmpi(get(aa(i),'type'),'uimenu')
        try
            set(aa(i),'BackgroundColor',c(i,:));
        end
    end
end
cd(wd)
