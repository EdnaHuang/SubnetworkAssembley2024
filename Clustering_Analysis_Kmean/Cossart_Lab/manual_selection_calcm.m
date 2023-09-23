function varargout = manual_selection_calcm(varargin)
% MANUAL_SELECTION_CALCM MATLAB code for manual_selection_calcm.fig
%      MANUAL_SELECTION_CALCM, by itself, creates a new MANUAL_SELECTION_CALCM or raises the existing
%      singleton*.
%
%      H = MANUAL_SELECTION_CALCM returns the handle to a new MANUAL_SELECTION_CALCM or the handle to
%      the existing singleton*.
%
%      MANUAL_SELECTION_CALCM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MANUAL_SELECTION_CALCM.M with the given input arguments.
%
%      MANUAL_SELECTION_CALCM('Property','Value',...) creates a new MANUAL_SELECTION_CALCM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before manual_selection_calcm_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to manual_selection_calcm_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help manual_selection_calcm

% Last Modified by GUIDE v2.5 21-Jul-2017 14:07:33

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @manual_selection_calcm_OpeningFcn, ...
                   'gui_OutputFcn',  @manual_selection_calcm_OutputFcn, ...
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


% --- Executes just before manual_selection_calcm is made visible.
function manual_selection_calcm_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to manual_selection_calcm (see VARARGIN)

% Choose default command line output for manual_selection_calcm
handles.output = hObject;
handles.k=[];
handles.xnew=[];
handles.ynew=[];
handles.x=[];
handles.y=[];
handles.FileName=[];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes manual_selection_calcm wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = manual_selection_calcm_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
a1=round(get(hObject,'Value'));
imshow(histeq(mat2gray(handles.k(:,:,a1))))
hold on;
scatter(handles.xnew,handles.ynew,'filled','b')
guidata(hObject, handles);




% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.FileName,handles.PathName,handles.FilterIndex]=uigetfile('*.mat');
addpath(handles.PathName);
%replace with Y if M2 undefined variable
load(handles.FileName,'M5');
handles.k=M5;
sz = size(handles.k);
imshow(histeq(mat2gray(handles.k(:,:,1))))
handles.xnew=[];
handles.ynew=[];
handles.x=[];
handles.y=[];
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold on;
[handles.x, handles.y]=ginput();
handles.xnew=[handles.xnew handles.x'];
handles.ynew=[handles.ynew handles.y'];
scatter(handles.xnew,handles.ynew,'filled','b')
guidata(hObject, handles);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
SSS='sel.mat';
sel=[];
sel(:,1)=handles.xnew;
sel(:,2)=handles.ynew;
save(SSS,'sel')
guidata(hObject, handles);


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold off;
imshow(imadjust(mat2gray(handles.k(:,:,1))))
guidata(hObject, handles);
