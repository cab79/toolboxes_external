function [varargout] = EP_den(varargin)
% EP_DEN M-file for EP_den.fig
%      EP_DEN, by itself, creates a new EP_DEN or raises the existing
%      singleton*.
%
%      H = EP_DEN returns the handle to a new EP_DEN or the handle to
%      the existing singleton*.
%
%      EP_DEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EP_DEN.M with the given input arguments.
%
%      EP_DEN('Property','Value',...) creates a new EP_DEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EP_den_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EP_den_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EP_den

% Last Modified by GUIDE v2.5 23-Nov-2012 14:32:09

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EP_den_OpeningFcn, ...
                   'gui_OutputFcn',  @EP_den_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State,varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before EP_den is made visible.
function EP_den_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EP_den (see VARARGIN)

% Choose default command line output for EP_den
handles.output = hObject;
handles.datatype ='ASCII';
set(handles.scales_popupmenu,'Value',5)
handles.par.samples = 512;
handles.par.stim = 129;
handles.par.sr = 512;
handles.par.scales = 5;
handles.par.plot_type = 'coeff';
handles.par.max_trials = 100;
handles.par.max_contour = 100;
handles.par.t_min=0.35;
handles.par.t_max=0.7;

set(handles.ncontour_edit,'value',40);
set(handles.nsingle_trials_edit,'value',10);
addpath(pwd);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes EP_den wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EP_den_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Load_data_button.
function Load_data_button_Callback(hObject, eventdata, handles)
% hObject    handle to Load_data_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.Denoise_Neigh_button,'value',0);
set(handles.Denoise_NZT_button,'value',0);
set(handles.Select_Scale_Coeffs_button,'value',0);
set(handles.load_den_coeff_button,'value',0);
set(handles.Plot_coefficients_button,'value',1);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Latency_Corrected_Average_button,'value',0);
handles.par.plot_type = 'coeff';

switch char(handles.datatype)
    case 'ASCII'  
        [filename, pathname] = uigetfile('*.asc','Select file');
        set(handles.file_name,'string',['Loading:    ' pathname filename]);
        eval(['cd ' pathname]);
        x=load([pathname filename]);                      % Load data
end
handles.par.filename = filename;
guidata(hObject, handles);

%Do average ERP.
samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;
x=x(:);
sweeps =length(x)/samples;
xx=reshape(x,samples,sweeps)';
av=mean(xx,1);

[coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);

axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av,'color',[0.5 0.5 0.5])
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])    
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

USER_DATA = get(handles.EP_auto_den_figure,'userdata');
USER_DATA{1} = [];   %denoised coefficients
USER_DATA{2} = x;    % single trials 
USER_DATA{3} = xx;  % reshaped single trials
USER_DATA{4} = coeff;  % wavelet decomposition coefficints        
USER_DATA{5} = [];   %denoised average EP (denav)
USER_DATA{6}=av;    % average EP
USER_DATA{7}=[]; % latency corrected average
USER_DATA{8}=[];     % YDEN
USER_DATA{9}=yo;   % original bands
USER_DATA{10}=[];   % denoised bands
USER_DATA{11}=[];  % current scale
%USER_DATA{12}=xmin; % time range for manual den
%USER_DATA{13}=xmax; % time range for manual den
USER_DATA{14}=[];   % single trials amplitudes (to remove the rectangle from the contour plot)
USER_DATA{15}=[];  % single trials latencies
%USER_DATA{16}=N_sc; %number of scales
%USER_DATA{16}=amp_av;
%USER_DATA{17}=lat_av;
USER_DATA{18}=[]; % latency corrected single trials
USER_DATA{19}=[];
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

%Does plot for the first time
axes(handles.axes2)
cla
set(handles.axes2,'xlimmode','auto');
set(handles.axes2,'ylimmode','auto');
step = 1/(sc+2):1/(sc+2):1;
for i=1:sc+1
    scaling_factor = 1.5 * max(max(abs(coeff(i,:)))) * (sc+1);
    aux = coeff(i,:)/ scaling_factor;
    plot(((1:samples)-stim+1)/sr,aux+step(sc+2-i),'color', [0.5 0.5 0.5])
    hold on
end
xlim([(1-stim+1)/sr (samples-stim+1)/sr])   
xrange = get(handles.axes2,'xlim') * 1.1;          
line([0 0],[0.05 0.95],'Linestyle',':','color','k')
for i=1:sc
    texto =['D' num2str(i)];
    text(xrange(1),step(sc+2-i)+0.01,texto);
end
texto =['A' num2str(sc)];
text(xrange(1),step(1)+0.01,texto);
axis off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Denoise_Neigh_button.
function Denoise_Neigh_button_Callback(hObject, eventdata, handles)
% hObject    handle to Denoise_Neigh_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Denoise_Neigh_button
set(handles.Denoise_Neigh_button,'value',1);
set(handles.Denoise_NZT_button,'value',0);
set(handles.Select_Scale_Coeffs_button,'value',0);
set(handles.load_den_coeff_button,'value',0);
set(handles.add_button,'value',0);
set(handles.Latency_Corrected_Average_button,'value',0);

USER_DATA = get(handles.EP_auto_den_figure,'userdata');
samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;
x=USER_DATA{2};
av=USER_DATA{6};

[coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);

USER_DATA{1}=den_coeff;
USER_DATA{5}=denav;
USER_DATA{10}=y;
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

handles.par.den_coeff=den_coeff;
guidata(hObject, handles);

%plot_average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

[YDEN]=st_den(x,den_coeff,handles);
USER_DATA{8}=YDEN;
USER_DATA{14}=[];   % single trials amplitudes (to remove the rectangle from the contour plot)
USER_DATA{18}=[];  % remove the latenct corrected single trials 
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

%plot coefficients
Plot_coefficients(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Denoise_NZT_button.
function Denoise_NZT_button_Callback(hObject, eventdata, handles)
% hObject    handle to Denoise_NZT_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Denoise_NZT_button
set(handles.Denoise_Neigh_button,'value',0);
set(handles.Denoise_NZT_button,'value',1);
set(handles.Select_Scale_Coeffs_button,'value',0);
set(handles.load_den_coeff_button,'value',0);
set(handles.add_button,'value',0);
set(handles.Latency_Corrected_Average_button,'value',0);

samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;
USER_DATA = get(handles.EP_auto_den_figure,'userdata');
x=USER_DATA{2};
av=USER_DATA{6};

[coeff,denav,den_coeff,y]= Run_NZT(av,handles);

USER_DATA{1}=den_coeff;
USER_DATA{5}=denav;
USER_DATA{10}=y;
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

handles.par.den_coeff=den_coeff;
guidata(hObject, handles);

%plot average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

[YDEN]=st_den(x,den_coeff,handles);
USER_DATA{8}=YDEN;
USER_DATA{14}=[];    % single trials amplitudes (to remove the rectangle from the contour plot)
USER_DATA{18}=[];  % remove the latency corrected single trials 
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

%plot coefficients
Plot_coefficients(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on selection change in scales_popupmenu.
function scales_popupmenu_Callback(hObject, eventdata, handles)
% hObject    handle to scales_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns scales_popupmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from
%        scales_popupmenu

sc_old = handles.par.scales;
sc = get(hObject,'Value');
if 2^sc <= handles.par.samples
    handles.par.scales = sc;
    guidata(hObject, handles);
else 
    errordlg('Too many scales for the Samples size');
    set(handles.scales_popupmenu,'Value',sc_old);
end

if strcmp(handles.par.plot_type,'single') ==1 | strcmp(handles.par.plot_type,'contour')==1;
    handles.par.plot_type = 'coeff';
    guidata(hObject, handles);
    set(handles.Plot_coefficients_button,'value',1);
    set(handles.Plot_bands_button,'value',0);
    set(handles.Plot_single_trials_button,'value',0);
    set(handles.Plot_contour_plot_button,'value',0);
end
set(handles.Denoise_Neigh_button,'value',0);
set(handles.Denoise_NZT_button,'value',0);

USER_DATA = get(handles.EP_auto_den_figure,'userdata');
USER_DATA{1} = [];
USER_DATA{8}=[];
USER_DATA{10}=[];   % denoised bands
if ~isempty(USER_DATA{3}) 
    sr = handles.par.sr;
    stim = handles.par.stim;
    samples = handles.par.samples;
    x = USER_DATA{2};                                           
    xx = USER_DATA{3};
    av = mean(xx,1);

    [coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);

    %plot_average
    axes(handles.average_plot)
    cla
    plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
    yrange = get(handles.average_plot,'ylim');
    line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
    xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
    title('Average ERP','fontsize',14)
    xlabel('Time (sec)','fontsize',10)

    %USER_DATA{1} = den_coeff;
    USER_DATA{2} = x;
    USER_DATA{3} = xx;
    USER_DATA{4} = coeff;          
    USER_DATA{6}=av;
    USER_DATA{9}=yo;
    set(handles.EP_auto_den_figure,'userdata',USER_DATA)
    
    %plot_coefficients
    axes(handles.axes2)
    cla
    set(handles.axes2,'xlimmode','auto');
    set(handles.axes2,'ylimmode','auto');
    step = 1/(sc+2):1/(sc+2):1;

    for i=1:sc+1
        scaling_factor1 = 1.5 * max(max(abs(coeff(i,:)))) * (sc+1);
        aux1 = coeff(i,:)/ scaling_factor1;
        plot(((1:samples)-stim+1)/sr,aux1+step(sc+2-i),'color', [0.5 0.5 0.5])
    end
    xlim([(1-stim+1)/sr (samples-stim+1)/sr])          
    xrange = get(handles.axes2,'xlim') * 1.1;          
    line([0 0],[0.05 0.95],'Linestyle',':','color','k')
    for i=1:sc
        texto =['D' num2str(i)];
        text(xrange(1),step(sc+2-i)+0.01,texto);
    end
    texto =['A' num2str(sc)];
    text(xrange(1),step(1)+0.01,texto);
    axis off
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function samples_edit_Callback(hObject, eventdata, handles)
% hObject    handle to samples_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of samples_edit as text
%        str2double(get(hObject,'String')) returns contents of samples_edit as a double
samples_old = handles.par.samples;
sc = handles.par.scales;
stim = handles.par.stim;
samples = str2num(get(handles.samples_edit,'string'));
if any(factor(samples)-2) == 0  
    if 2^sc <= samples
        if stim <= samples 
            handles.par.samples = samples;
            guidata(hObject, handles);
        else
            errordlg('Correct Stim. first (should be lower than Samples)');
            set(handles.samples_edit,'Value',samples_old);
            set(handles.samples_edit,'String',num2str(samples_old));
        end
    else 
        errordlg('Too few samples for the number of scales');
        set(handles.samples_edit,'Value',samples_old);
        set(handles.samples_edit,'String',num2str(samples_old));
    end
else 
    errordlg('Number of samples should be a power of 2');
    set(handles.samples_edit,'Value',samples_old);
    set(handles.samples_edit,'String',num2str(samples_old));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stim_edit_Callback(hObject, eventdata, handles)
% hObject    handle to stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of stim_edit as text
%        str2double(get(hObject,'String')) returns contents of stim_edit as a double
stim_old = handles.par.stim;
samples = handles.par.samples;
stim = str2num(get(handles.stim_edit,'string'));
if stim <= samples & stim > 0
    handles.par.stim = stim;
    guidata(hObject, handles);
else
    errordlg('Stim. should be lower than Samples');
    set(handles.stim_edit,'Value',stim_old);
    set(handles.stim_edit,'String',num2str(stim_old));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sr_edit_Callback(hObject, eventdata, handles)
% hObject    handle to sr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of sr_edit as text
%        str2double(get(hObject,'String')) returns contents of sr_edit as a double

sr=str2num(get(handles.sr_edit,'string'));
handles.par.sr=sr;
guidata(hObject,handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Save_button.
function Save_button_Callback(hObject, eventdata, handles)
% hObject    handle to Save_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

USER_DATA=get(handles.EP_auto_den_figure,'userdata');
 den_coeff=USER_DATA{1};
 x=USER_DATA{2};
 xx=USER_DATA{3};
 coeff=USER_DATA{4};          
 denav=USER_DATA{5};
 av_lc=USER_DATA{7};
 YDEN=USER_DATA{8};
%  yo=USER_DATA{9};
%  y=USER_DATA{10};
% tmin=USER_DATA{12};
% tmax=USER_DATA{13};
amp=USER_DATA{14};
lat=USER_DATA{15};
% amp_av=USER_DATA{16};
% lat_av=USER_DATA{17};
xx_lc=USER_DATA{18};
t_min=handles.par.t_min;
t_max=handles.par.t_max;
outfile=[handles.par.filename '_den.mat'];
save(outfile,'den_coeff','x','xx','coeff','denav','av_lc','YDEN','amp','lat','xx_lc','t_min','t_max');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Plot_coefficients_button.
function Plot_coefficients_button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_coefficients_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% Plot_coefficients_button

set(handles.Plot_coefficients_button,'value',1);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',0);
handles.par.plot_type = 'coeff';
guidata(hObject, handles);

Plot_coefficients(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Plot_contour_plot_button.
function Plot_contour_plot_button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_contour_plot_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% Plot_contour_plot_button

set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',1);
handles.par.plot_type = 'contour';
guidata(hObject, handles);

Plot_coefficients(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in Plot_bands_button.
function Plot_bands_button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_bands_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Plot_bands_button

set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',1);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',0);
handles.par.plot_type = 'bands';
guidata(hObject, handles);

Plot_coefficients(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Plot_single_trials_button.
function Plot_single_trials_button_Callback(hObject, eventdata, handles)
% hObject    handle to Plot_single_trials_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% Plot_single_trials_button

set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',1);
set(handles.Plot_contour_plot_button,'value',0);
handles.par.plot_type = 'single';
guidata(hObject, handles);

Plot_coefficients(handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ncontour_edit_Callback(hObject, eventdata, handles)
% hObject    handle to ncontour_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ncontour_edit as text
%        str2double(get(hObject,'String')) returns contents of ncontour_edit as a double
max_contour = str2num(get(handles.ncontour_edit,'string'));
handles.par.max_contour = max_contour;
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nsingle_trials_edit_Callback(hObject, eventdata, handles)
% hObject    handle to nsingle_trials_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of nsingle_trials_edit as text
%        str2double(get(hObject,'String')) returns contents of
%        nsingle_trials_edit as a double

max_trials = str2num(get(handles.nsingle_trials_edit,'string'));
handles.par.max_trials = max_trials;
guidata(hObject, handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in Select_Scale_Coeffs_button.
function Select_Scale_Coeffs_button_Callback(hObject, eventdata, handles)
% hObject    handle to Select_Scale_Coeffs_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Select_Scale_Coeffs_button
set(handles.Denoise_Neigh_button,'value',0);
set(handles.Denoise_NZT_button,'value',0);
set(handles.Select_Scale_Coeffs_button,'value',1);
set(handles.load_den_coeff_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Plot_coefficients_button,'value',1);
set(handles.Latency_Corrected_Average_button,'value',0);
handles.par.plot_type = 'coeff';
guidata(hObject, handles);               
Plot_coefficients(handles);
sc=handles.par.scales;
sr = handles.par.sr;
stim=handles.par.stim;
samples=handles.par.samples;
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
denav=USER_DATA{5};
av=USER_DATA{6};
c_sc=USER_DATA{11};


step = 1/(sc+2):1/(sc+2):1;
axes(handles.axes2)
rect=getrect;
ymin=rect(1,2);
ymax=ymin+rect(1,4);
for i=1:sc+1
     if (step(i)<=ymax) && (ymax<=step(i+1))
         c_sc=sc+2-i;                                     % current scale
     end
end

xmin=rect(1,1);                     
xmax=xmin+rect(1,3);
xmin=round(xmin*sr+stim);
xmax=round(xmax*sr+stim);
if xmin<0
    xmin=1;  
end
if xmax>samples
    xmax=samples;
end

N_sc=0;                                 % define the number of selected scales
for i=1:length(step)
if (ymin<step(i) )&&( step(i)<ymax)
N_sc=N_sc+1;
end
end

axes(handles.average_plot)       
% rh=findobj('type','rectangle');  % remove rectangle in case of having find peaks rect
% delete(rh)
%plot_average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

USER_DATA{11}=c_sc;
USER_DATA{12}=xmin;
USER_DATA{13}=xmax;
USER_DATA{16}=N_sc;
USER_DATA{14}=[];  % single trials amplitudes (to remove the rectangle from the contour plot)
USER_DATA{18}=[];  % remove the latenct corrected single trials 
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in add_button.
function add_button_Callback(hObject, eventdata, handles)
% hObject    handle to add_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of add_button
set(handles.Select_Scale_Coeffs_button,'value',1);
set(handles.add_button,'value',1);
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
den_coeff=USER_DATA{1};
x=USER_DATA{2};
coeff=USER_DATA{4}; 
av=USER_DATA{6};
c_sc=USER_DATA{11};
xmin=USER_DATA{12};
xmax=USER_DATA{13};
N_sc=USER_DATA{16};

samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;
if isempty(den_coeff)
        den_coeff=zeros(sc+1,samples);
end

for j=1:N_sc    
    for i=xmin:xmax
            den_coeff(c_sc+(j-1),i)=coeff(c_sc+(j-1),i);
    end   
end
USER_DATA{1}=den_coeff;

[denav,y]=st_den(av,den_coeff,handles);
USER_DATA{5}=denav;
USER_DATA{10}=y;

[YDEN]=st_den(x,den_coeff,handles);
USER_DATA{8}=YDEN;
USER_DATA{11}=[];
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

handles.par.den_coeff=den_coeff;
guidata(hObject, handles);

Plot_coefficients(handles);
%plot average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

set(handles.add_button,'value',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in remove_button.
function remove_button_Callback(hObject, eventdata, handles)
% hObject    handle to remove_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of remove_button
set(handles.Select_Scale_Coeffs_button,'value',1);
set(handles.remove_button,'value',1);
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
den_coeff=USER_DATA{1};
x=USER_DATA{2};
coeff=USER_DATA{4}; 
av=USER_DATA{6};
c_sc=USER_DATA{11};
xmin=USER_DATA{12};
xmax=USER_DATA{13};
N_sc=USER_DATA{16};

samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;

if isempty(den_coeff)
    den_coeff=zeros(sc+1,samples);
end

for j=1:N_sc  
for i=xmin:xmax
        den_coeff(c_sc+(j-1),i)=0;
end
end
USER_DATA{1}=den_coeff;


[denav,y]=st_den(av,den_coeff,handles);
USER_DATA{5}=denav;
USER_DATA{10}=y;

[YDEN]=st_den(x,den_coeff,handles);
USER_DATA{8}=YDEN;
USER_DATA{11}=[];
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

handles.par.den_coeff=den_coeff;
guidata(hObject, handles);

Plot_coefficients(handles);
%plot average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

set(handles.remove_button,'value',0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in load_den_coeff_button.
function load_den_coeff_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_den_coeff_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of load_den_coeff_button
set(handles.Denoise_Neigh_button,'value',0);
set(handles.Denoise_NZT_button,'value',0);
set(handles.Select_Scale_Coeffs_button,'value',0);
set(handles.load_den_coeff_button,'value',1);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',0);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Plot_coefficients_button,'value',1);
set(handles.Latency_Corrected_Average_button,'value',0);
handles.par.plot_type = 'coeff';

USER_DATA=get(handles.EP_auto_den_figure,'userdata');
x=USER_DATA{2};
coeff=USER_DATA{4};
av=USER_DATA{6};
%den_coeff=handles.par.den_coeff;
[filename, pathname] = uigetfile('*.mat','Select file');
matfile=load([pathname filename]);  
den_coeff=matfile.den_coeff;

samples=handles.par.samples;
sr=handles.par.sr;
stim=handles.par.stim;
sc=handles.par.scales;

[denav,y,den_coeff]=st_den(av,den_coeff,handles);
USER_DATA{5}=denav;
USER_DATA{10}=y;

[YDEN]=st_den(x,den_coeff,handles);
USER_DATA{8}=YDEN;
USER_DATA{1}=den_coeff;
set(handles.EP_auto_den_figure,'userdata',USER_DATA)

Plot_coefficients(handles);
%plot average
axes(handles.average_plot)
cla
plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
hold on
plot(((1:samples)-stim+1)/sr,denav, 'color','r')
yrange = get(handles.average_plot,'ylim');
line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
title('Average ERP','fontsize',14)
xlabel('Time (sec)','fontsize',10)

set(handles.add_button,'value',0);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function t_min_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_min_edit as text
%        str2double(get(hObject,'String')) returns contents of t_min_edit as a double;
t_min_old=handles.par.t_min;
t_min=str2num(get(handles.t_min_edit,'string'));
handles.par.plot_type = 'single';
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
YDEN=USER_DATA{8};

if isempty(YDEN) 
    errordlg('First do the denoising')
    set(handles.t_min_edit,'Value',t_min_old);
    set(handles.t_min_edit,'String',num2str(t_min_old));
else
    set(handles.Plot_coefficients_button,'value',0);
    set(handles.Plot_bands_button,'value',0);
    set(handles.Plot_single_trials_button,'value',1);
    set(handles.Plot_contour_plot_button,'value',0);
    handles.par.t_min = t_min;
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function t_max_edit_Callback(hObject, eventdata, handles)
% hObject    handle to t_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of t_max_edit as text
%        str2double(get(hObject,'String')) returns contents of t_max_edit as a double
t_max_old=handles.par.t_max;
t_max=str2num(get(handles.t_max_edit,'string'));
handles.par.plot_type = 'single';
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
denav=USER_DATA{5};
av=USER_DATA{6};
YDEN=USER_DATA{8};
t_min=handles.par.t_min;
sr=handles.par.sr;
stim=handles.par.stim;
samples=handles.par.samples;

if isempty(YDEN) 
    errordlg('First do the denoising')
    set(handles.t_max_edit,'Value',t_max_old);
    set(handles.t_max_edit,'String',num2str(t_max_old));
else
    set(handles.Plot_coefficients_button,'value',0);
    set(handles.Plot_bands_button,'value',0);
    set(handles.Plot_single_trials_button,'value',1);
    set(handles.Plot_contour_plot_button,'value',0);
    handles.par.t_max = t_max;
    guidata(hObject, handles);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Positive_radiobutton.
function Positive_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to Positive_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Positive_radiobutton
set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',1);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Positive_radiobutton,'value',1);
set(handles.Negative_radiobutton,'value',0);
handles.par.plot_type = 'single';
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
denav=USER_DATA{5};
av=USER_DATA{6};
YDEN=USER_DATA{8};
sr=handles.par.sr;
stim=handles.par.stim;
samples=handles.par.samples;
t_min=handles.par.t_min;
t_max=handles.par.t_max;

if isempty(YDEN) 
    errordlg('First do the denoising')
else
    tmin = ceil(t_min * sr + stim-1);
    USER_DATA{12}=tmin;
    tmax = floor(t_max * sr + stim-1);
    USER_DATA{13}=tmax;
    
    [amp_av lat_av]=max(denav(tmin:tmax));
    lat_av=(lat_av + tmin - stim) * 1000/ sr;
    amp = zeros(size(YDEN,1),1);
     lat = zeros(size(YDEN,1),1);
     for i = 1:size(YDEN,1)
            [amp(i) lat(i)] = max(YDEN(i,tmin:tmax));            
            lat(i) = (lat(i) + tmin - stim) * 1000/ sr;    %in ms
     end
     USER_DATA{14}=amp;
     USER_DATA{15}=lat;
     USER_DATA{16}=amp_av;
     USER_DATA{17}=lat_av;
     USER_DATA{18}=[];  % remove the latenct corrected single trials 
     % plot average
     lat_av=lat_av/1000;
    axes(handles.average_plot)
    cla
    hold on
    plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
    plot(((1:samples)-stim+1)/sr,denav, 'color','r')
    plot(lat_av,amp_av,'b:*')
    rectangle('position',[t_min,-max(av),t_max-t_min,2*max(av)],'edgecolor','b')
    yrange = get(handles.average_plot,'ylim');
    line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
    xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
    title('Average ERP','fontsize',14)
    xlabel('Time (sec)','fontsize',10)
set(handles.EP_auto_den_figure,'userdata',USER_DATA)
Plot_coefficients(handles);
end
%waitsecs(1);
set(handles.Positive_radiobutton,'value',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Negative_radiobutton.
function Negative_radiobutton_Callback(hObject, eventdata, handles)
% hObject    handle to Negative_radiobutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Negative_radiobutton
set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',1);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Positive_radiobutton,'value',0);
set(handles.Negative_radiobutton,'value',1);
handles.par.plot_type = 'single';
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
denav=USER_DATA{5};
av=USER_DATA{6};
YDEN=USER_DATA{8};
sr=handles.par.sr;
stim=handles.par.stim;
samples=handles.par.samples;
t_min=handles.par.t_min;
t_max=handles.par.t_max;

if isempty(YDEN) 
       errordlg('First do the denoising')
else

    tmin = ceil(t_min * sr + stim-1);
    USER_DATA{12}=tmin;
    tmax = floor(t_max * sr + stim-1);

    USER_DATA{13}=tmax;
    
    [amp_av lat_av]=min(denav(tmin:tmax));
    lat_av=(lat_av + tmin - stim) * 1000/ sr;
    amp = zeros(size(YDEN,1),1);
    lat = zeros(size(YDEN,1),1);
     for i = 1:size(YDEN,1)
            [amp(i) lat(i)] = min(YDEN(i,tmin:tmax));            
            lat(i) = (lat(i) + tmin - stim) * 1000/ sr;    %in ms
     end
      USER_DATA{14}=amp;
      USER_DATA{15}=lat;
      USER_DATA{16}=amp_av;
      USER_DATA{17}=lat_av;
      USER_DATA{18}=[];  % remove the latency corrected single trials 
      % plot average
    lat_av=lat_av/1000;
    axes(handles.average_plot)
    cla
    hold on
    plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
    plot(((1:samples)-stim+1)/sr,denav, 'color','r')
    plot(lat_av,amp_av,'b:*')
    rectangle('position',[t_min,-max(av),t_max-t_min,2*max(av)],'edgecolor','b')
    yrange = get(handles.average_plot,'ylim');
    line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
    xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
    title('Average ERP','fontsize',14)
    xlabel('Time (sec)','fontsize',10)
set(handles.EP_auto_den_figure,'userdata',USER_DATA)
Plot_coefficients(handles);
end
%waitsecs(0.5);
set(handles.Negative_radiobutton,'value',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes on button press in Latency_Corrected_Average_button.
function Latency_Corrected_Average_button_Callback(hObject, eventdata, handles)
% hObject    handle to Latency_Corrected_Average_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of
% Latency_Corrected_Average_button
set(handles.Plot_coefficients_button,'value',0);
set(handles.Plot_bands_button,'value',0);
set(handles.Plot_single_trials_button,'value',1);
set(handles.Plot_contour_plot_button,'value',0);
set(handles.Latency_Corrected_Average_button,'value',1);
handles.par.plot_type = 'single';
USER_DATA=get(handles.EP_auto_den_figure,'userdata');
x=USER_DATA{2};
denav=USER_DATA{5};
av=USER_DATA{6};
YDEN=USER_DATA{8};
amp=USER_DATA{14};
lat=USER_DATA{15};
amp_av=USER_DATA{16};
lat_av=USER_DATA{17};
sr=handles.par.sr;
stim=handles.par.stim;
samples=handles.par.samples;
sweeps =length(x)/samples;
xx=reshape(x,samples,sweeps)';


if isempty(YDEN) || isempty(amp)
       errordlg('select the peak')
else
    for j=1:size(YDEN,1)
        n=round((lat(j)-lat_av)*sr/1000);
        lc=xx(j,:);   
        y=zeros(1,length(lc));
        for i=1:length(lc)
            if i+n >= 1 & i+n <= length(lc)
                y(i) = lc(i+n);
            end
        end   
        xx_lc(j,:)=y;
    end
    av_lc=mean(xx_lc);
    USER_DATA{18}=xx_lc;
    USER_DATA{7}=av_lc;
    
    %plot average
    axes(handles.average_plot)
    cla
    hold on
    plot(((1:samples)-stim+1)/sr,av, 'color',[0.5 0.5 0.5])
     plot(((1:samples)-stim+1)/sr,denav, 'color','r')
    plot(((1:samples)-stim+1)/sr,av_lc, 'color','b')
    yrange = get(handles.average_plot,'ylim');
    line([0 0],[yrange(1) yrange(2)],'Linestyle',':','color','k')
    xlim([(1-stim+1)/sr (samples-stim+1)/sr])         
    title('Average ERP','fontsize',14)
    xlabel('Time (sec)','fontsize',10)
    
     av=av_lc;
    [coeff,denav,den_coeff,y,yo]= Run_Neigh(av,handles);
    yo_lc=yo;
    USER_DATA{19}=yo_lc;
    set(handles.EP_auto_den_figure,'userdata',USER_DATA)
    Plot_coefficients(handles);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --- Executes during object creation, after setting all properties.
function samples_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to samples_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function scales_popupmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to scales_popupmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function sr_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sr_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function stim_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stim_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function nsingle_trials_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to nsingle_trials_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function ncontour_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ncontour_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function data_type_popmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to data_type_popmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function t_min_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_min_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes during object creation, after setting all properties.
function t_max_edit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to t_max_edit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










