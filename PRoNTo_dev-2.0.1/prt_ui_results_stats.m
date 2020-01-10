function varargout = prt_ui_results_stats(varargin)
% PRT_UI_RESULTS_STATS MATLAB code for prt_ui_results_stats.fig
%
% PRT_UI_RESULTS_STATS, by itself, creates a new PRT_UI_RESULTS_STATS or raises the
% existing singleton*.
%
% H = PRT_UI_RESULTS_STATS returns the handle to a new PRT_UI_RESULTS_STATS or the
% handle to the existing singleton*.
%
% PRT_UI_RESULTS_STATS('CALLBACK',hObject,eventData,handles,...) calls the local
% function named CALLBACK in PRT_UI_RESULTS_STATS.M with the given input arguments.
%
% PRT_UI_RESULTS_STATS('Property','Value',...) creates a new PRT_UI_RESULTS_STATS or
% raises the existing singleton*.  Starting from the left, property value
% pairs are applied to the GUI before prt_ui_results_stats_OpeningFcn gets called.
% An unrecognized property name or invalid value makes property application
% stop.  All inputs are passed to prt_ui_results_stats_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id: $

% Edit the above text to modify the response to help prt_ui_results_stats

% Last Modified by GUIDE v2.5 26-Jan-2015 17:01:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @prt_ui_results_stats_OpeningFcn, ...
    'gui_OutputFcn',  @prt_ui_results_stats_OutputFcn, ...
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


% --- Executes just before prt_ui_results_stats is made visible.
function prt_ui_results_stats_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_results_stats (see VARARGIN)

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
set(handles.save_perm_weights,'Visible','off')
% Choose default command line output for prt_ui_results_stats
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_results_stats wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_results_stats_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

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
statsbutton_Callback(hObject, eventdata, handles);

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

nplot         = get(handles.plotmenu,'String');
plotm         = get(handles.plotmenu,'Value');
plotchosen    = num2str(plotm);
if plotm>length(nplot)  % reset to 1 if list of available plot smaller than chosen plot
    set(handles.plotmenu,'Value',1);
    plotm = 1;
    plotchosen    = num2str(plotm);
end
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
pos = [0.108 0.1101 0.8543 0.8165];
set(handles.axes5,'Position',pos);


if strcmp(PRT.model(model).input.type,'classification')
    %     mclass        = length(handles.PRT.model(model).output.stats.c_acc);
    %     if mclass ~= 2
    %         multiplot  = 4;
    %         plotchosen = num2str(multiplot(plotm));
    %     end
    
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

% % Check if stats exist
% % -------------------------------------------------------------------------
% if isfield(handles, 'stats')
%     stats = handles.stats;
% end

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
    
end

%%% SHOW STATS
    
switch stats.type
       
    case 'class'
        
        % Disable regression
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
        
        set(handles.accuracytext,'String','Total accuracy:','Visible','on');
        set(handles.baccuracytext,'String','Balanced accuracy (BA):','Visible','on');
        set(handles.classaccuracytext,'String','Class accuracy (CA):','Visible','on');
        set(handles.ppvtext,'String','Class predictive value: ','Visible','on');
        
        set(handles.acctext,'String',sprintf('%3.2f %%',stats.macc*100),'Visible','on');
        set(handles.bacctext,'String',sprintf('%3.2f %%',stats.mbacc*100),'Visible','on');
        set(handles.cacctext,'String',sprintf('  %3.2f %%',stats.mcacc*100),'Visible','on');
        set(handles.cpvval,'String',sprintf('  %3.2f %%',...
            stats.mcpv*100),'Visible','on');
        
        set(handles.pvalbacc,'String','BA p-value:','Visible','on');
        set(handles.pvalcacc,'String','CA p-value:','Visible','on');
        set(handles.pbacc,'String','N. A.','Visible','on');
        set(handles.pcacc,'String','N. A.','Visible','on');
        
        if isfield(stats,'show_perm')
            
            if stats.show_perm

                set(handles.pbacc,'String',sprintf(' %3.4f',stats.perm.pvalue_b_acc));
                set(handles.pcacc,'String',sprintf(' %3.4f',stats.perm.pvalue_c_acc));
                
            end
            
        end
         
    case 'reg'
        
        set(handles.accuracytext,'Visible','off');
        set(handles.baccuracytext,'Visible','off');
        set(handles.classaccuracytext,'Visible','off');
        set(handles.ppvtext,'Visible','off');
        set(handles.acctext,'Visible','off');
        set(handles.bacctext,'Visible','off');
        set(handles.cacctext,'Visible','off');
        set(handles.cpvval,'Visible','off');
        set(handles.pvalbacc,'Visible','off');
        set(handles.pvalcacc,'Visible','off');
        set(handles.pbacc,'Visible','off');
        set(handles.pcacc,'Visible','off');
        
        set(handles.corrtext,'String','Correlation:','Visible','on');
        set(handles.corrvaltext,'String',sprintf('%3.2f',stats.corr),'Visible','on');
        
        if isfield(stats,'r2')
            set(handles.r2text,'String','Coeff. of determination (R2):','Visible','on');
            set(handles.r2valtext,'String',sprintf('%3.2f',stats.r2),'Visible','on');
            set(handles.pr2,'Visible','on','String','R2 p-value: N. A.');
        end
        
        set(handles.msetext,'String','MSE:','Visible','on');
        set(handles.msevaltext,'String',sprintf('%3.2f',stats.mse),'Visible','on');
        
        if isfield(stats,'nmse')
            set(handles.nmsetext,'String','Norm. MSE:','Visible','on');
            set(handles.nmsevaltext,'String',sprintf('%3.2f',stats.nmse),'Visible','on');
            set(handles.pnmse,'Visible','on','String','Norm. MSE p-value: N. A.');
        end
        
        set(handles.pcorr,'Visible','on','String','Correlation p-value: N. A.');   
        set(handles.pmse,'Visible','on','String','MSE p-value: N. A.');

        
        if isfield(stats,'show_perm')
            
            if stats.show_perm
                
                set(handles.pcorr,'Visible','on','String',sprintf('Correlation p-value: %3.4f',stats.perm.pval_corr));
                if isfield(stats.perm,'pval_r2'), set(handles.pr2,'Visible','on','String',sprintf('R2 p-value: %3.4f',stats.perm.pval_r2)); end
                set(handles.pmse,'Visible','on','String',sprintf('MSE p-value: %3.4f',stats.perm.pval_mse));
                if isfield(stats.perm,'pval_nmse'), set(handles.pnmse,'Visible','on','String',sprintf('Norm. MSE p-value: %3.4f',stats.perm.pval_nmse));end
                
            end
        end
end
    

stats.show_perm = 0;
handles.stats = stats;
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


% --- Executes on button press in showweights.
function showweights_Callback(hObject, eventdata, handles)
% hObject    handle to showweights (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in editPlotButton.
function editPlotButton_Callback(hObject, eventdata, handles)
% hObject    handle to editPlotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get stuff from the displayed figure
old_fig = gcf;
figure_children = get(old_fig,'Children');
children_axes = findall(figure_children,'Type','axes'); %important to get the legends too

% Create a new figure
fig_out = figure;
if length(children_axes) > 1 % There are legends in the figure
    axes_out = copyobj([children_axes(1); children_axes(2)], fig_out);
else % There are no legends in the figure
    axes_out = copyobj(children_axes, fig_out);
end
