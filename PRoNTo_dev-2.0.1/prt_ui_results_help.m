function varargout = prt_ui_results_help(varargin)
% PRT_UI_RESULTS_HELP MATLAB code for prt_ui_results_help.fig
% 
% PRT_UI_RESULTS_HELP, by itself, creates a new PRT_UI_RESULTS_HELP or 
% raises the existing singleton*.
%
% H = PRT_UI_RESULTS_HELP returns the handle to a new PRT_UI_RESULTS_HELP 
% or the handle to the existing singleton*.
%
% PRT_UI_RESULTS_HELP('CALLBACK',hObject,eventData,handles,...) calls the 
% local function named CALLBACK in PRT_UI_RESULTS_HELP.M with the given 
% input arguments.
%
% PRT_UI_RESULTS_HELP('Property','Value',...) creates a new 
% PRT_UI_RESULTS_HELP or raises the existing singleton*.  Starting from the
% left, property value pairs are applied to the GUI before 
% prt_ui_results_help_OpeningFcn gets called.  An unrecognized property 
% name or invalid value makes property application stop.  All inputs are 
% passed to prt_ui_results_help_OpeningFcn via varargin.
%
% *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
% instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%__________________________________________________________________________
% Copyright (C) 2011 Machine Learning & Neuroimaging Laboratory

% Written by M. J. Rosa
% $Id$

% Edit the above text to modify the response to help prt_ui_results_help

% Last Modified by GUIDE v2.5 13-Mar-2012 10:41:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @prt_ui_results_help_OpeningFcn, ...
                   'gui_OutputFcn',  @prt_ui_results_help_OutputFcn, ...
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


% --- Executes just before prt_ui_results_help is made visible.
function prt_ui_results_help_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to prt_ui_results_help (see VARARGIN)

% Choose default command line output for prt_ui_results_help
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes prt_ui_results_help wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = prt_ui_results_help_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
