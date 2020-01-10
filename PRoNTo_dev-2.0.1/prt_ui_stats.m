function varargout = prt_ui_stats(varargin)
% PRT_UI_STATS MATLAB code for prt_ui_stats.fig
%
% PRT_UI_STATS, by itself, creates a new PRT_UI_STATS or raises the
% existing singleton*.
%
% H = PRT_UI_STATS returns the handle to a new PRT_UI_STATS or the handle
% to the existing singleton*.
%
% PRT_UI_STATS('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_STATS.M with the given input arguments.
%
% PRT_UI_STATS('Property','Value',...) creates a new PRT_UI_STATS or raises
% the existing singleton*.  Starting from the left, property value pairs
% are applied to the GUI before prt_ui_stats_OpeningFcn gets called.  An
% unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_stats_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id$

% Edit the above text to modify the response to help prt_ui_stats

% Last Modified by GUIDE v2.5 29-Mar-2012 11:29:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prt_ui_stats_OpeningFcn, ...
    'gui_OutputFcn',  @prt_ui_stats_OutputFcn, ...
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


% --- Executes just before prt_ui_stats is made visible.
function prt_ui_stats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_stats (see VARARGIN)

%if window already exists, just put it as the current figure
Tag='stats';
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
    % set(handles.figure1,'Name','PRoNTo :: Stats table)
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
    
    rotate3d off
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
                        set(bb(j),'BackgroundColor',color.bg2)
                        cc=get(bb(j),'children');
                        for k=1:length(cc)
                            if strcmpi(get(cc(k),'type'),'uicontrol') && ...
                                    ~isempty(find(strcmpi(get(cc(k),'Style'),{'text',...
                                    'radiobutton','checkbox'})))
                                set(cc(k),'BackgroundColor',color.bg2)
                            elseif strcmpi(get(cc(k),'type'),'uicontrol') && ...
                                    ~isempty(find(strcmpi(get(cc(k),'Style'),'pushbutton')))
                                set(cc(k),'BackgroundColor',color.fr)
                            end
                            set(cc(k),'FontUnits','pixel')
                            xf=get(cc(k),'FontSize');
                            set(cc(k),'FontSize',ceil(FS*xf),'FontName',PF,...
                                'FontUnits','normalized','Units','normalized')
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
end

if ~isempty(varargin)
    
    stats = varargin{1};
    
    switch stats.type
        
        
        case 'class'
            
            set(handles.corrtext,'Visible','off');
            set(handles.corrvaltext,'Visible','off');
            
            set(handles.r2text,'Visible','off');
            set(handles.r2valtext,'Visible','off');
            
            set(handles.msetext,'Visible','off');
            set(handles.msevaltext,'Visible','off');
            
            set(handles.nmsetext,'Visible','off');
            set(handles.nmsevaltext,'Visible','off');
            
            set(handles.pcorr,'Visible','off');
            set(handles.pr2,'Visible','off');
            set(handles.pmse,'Visible','off');
            set(handles.pnmse,'Visible','off');
            
            set(handles.accuracytext,'String','Accuracy (acc):','Visible','on');
            set(handles.baccuracytext,'String','Balanced acc:','Visible','on');
            set(handles.classaccuracytext,'String','Class acc (%):','Visible','on');
            set(handles.ppvtext,'String','Class pv (%): ','Visible','on');
            set(handles.cpvval,'String',sprintf(' %3.1f ',...
                stats.mcpv*100),'Visible','on');
            set(handles.npvtext,'String','NPV: not yet available','Visible','off');
            
            set(handles.acctext,'String',sprintf('%3.1f %%',stats.macc*100),'Visible','on');
            set(handles.bacctext,'String',sprintf('%3.1f %%',stats.mbacc*100),'Visible','on');
            
            set(handles.cacctext,'String',sprintf(' %3.1f',stats.mcacc*100),'Visible','on');
            
            set(handles.pcorr, 'Visible','off');
            set(handles.pr2, 'Visible','off');
            set(handles.pmse,'Visible','off');
            set(handles.pnmse,'Visible','off');
            set(handles.pbacc,'Visible','off');
            set(handles.pcacc,'Visible','off');
            
            
            if isfield(stats,'show_perm')
                
                if stats.show_perm
                    
                    beep;
                    disp('...')
                    disp('Permutations results:')
                    disp(sprintf('Balanced accuracy p-value: %3.4f',stats.perm.pvalue_b_acc));
                    disp(sprintf('Class accuracy p-value:'));
                    disp(sprintf(' %3.4f',stats.perm.pvalue_c_acc));
                    
                end
                
            end
            
        case 'reg'
            
            set(handles.accuracytext,'String','Accuracy (acc):','Visible','off');
            set(handles.baccuracytext,'String','Balanced acc:','Visible','off');
            set(handles.classaccuracytext,'String','Class acc:','Visible','off');
            set(handles.pbacc, 'Visible','off');
            set(handles.pcacc,'Visible','off');
            
            set(handles.acctext,'Visible','off');
            set(handles.bacctext,'Visible','off');
            set(handles.cacctext,'Visible','off');
            
            set(handles.ppvtext,'Visible','off');
            set(handles.npvtext,'Visible','off');
            
            set(handles.corrtext,'String','Correlation:','Visible','on');
            set(handles.corrvaltext,'String',sprintf('%3.2f',stats.corr),'Visible','on');
            
            if isfield(stats,'r2')
                set(handles.r2text,'String','Coefficient of determination:','Visible','on');
                set(handles.r2valtext,'String',sprintf('%3.2f',stats.r2),'Visible','on');
            end
            
            set(handles.msetext,'String','MSE:','Visible','on');
            set(handles.msevaltext,'String',sprintf('%3.2f',stats.mse),'Visible','on');
            
            if isfield(stats,'nmse')
                set(handles.nmsetext,'String','Normalised MSE:','Visible','on');
                set(handles.nmsevaltext,'String',sprintf('%3.2f',stats.nmse),'Visible','on');
            end
            
            set(handles.pbacc, 'Visible','off');
            set(handles.pcacc,'Visible','off');
            set(handles.pcorr,'Visible','off');
            set(handles.pr2,'Visible','off');
            set(handles.pmse,'Visible','off');
            set(handles.pnmse,'Visible','off');
            
            if isfield(stats,'show_perm')
                
                if stats.show_perm
                    
                    beep;
                    disp('...')
                    disp('Permutations results:')
                    disp(sprintf('Correlation p-value: %3.4f',stats.perm.pval_corr));
                    if isfield(stats.perm,'pval_r2'), disp(sprintf('Coefficient of determination p-value: %3.4f',stats.perm.pval_r2)); end
                    disp(sprintf('Mean squared-error p-value: %3.4f',stats.perm.pval_mse));
                    if isfield(stats.perm,'pval_nmse'), disp(sprintf('Normalised mean squared-error p-value: %3.4f',stats.perm.pval_nmse)); end
                    
                end
            end
    end
    handles.prtdir=varargin{2};
    
end

rotate3d off
% Choose default command line output for prt_ui_stats
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_stats wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_stats_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;




% --------------------------------------------------------------------
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
    'Save figure as','Stats_table.png');
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
