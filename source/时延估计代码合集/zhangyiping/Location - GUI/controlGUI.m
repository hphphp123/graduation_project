function varargout = controlGUI(varargin)
% CONTROLGUI MATLAB code for controlGUI.fig
%      CONTROLGUI, by itself, creates a new CONTROLGUI or raises the existing
%      singleton*.
%
%      H = CONTROLGUI returns the handle to a new CONTROLGUI or the handle to
%      the existing singleton*.
%
%      CONTROLGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONTROLGUI.M with the given input arguments.
%
%      CONTROLGUI('Property','Value',...) creates a new CONTROLGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before controlGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to controlGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help controlGUI

% Last Modified by GUIDE v2.5 22-Jan-2016 13:23:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @controlGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @controlGUI_OutputFcn, ...
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


% --- Executes just before controlGUI is made visible.
function controlGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to controlGUI (see VARARGIN)

% Choose default command line output for controlGUI
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes controlGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = controlGUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
global Gtrigger Ttrigger
Gtrigger = 1;
Ttrigger = 0;

theta = (0:0.05:2)*pi;
x= cos(theta);
y = sin(theta);
LL = length(theta);
z = ones(LL);
plot3(x,y,z,'-b');
axis([-1.2 1.2 -1.2 1.2 0 2]);
xh = [-1.2  1.2];
yh = [0 0];
zh = [1 1];
hold on
plot3(xh,yh,zh,'-.b');
xv = [0 0];
yv = [-1.2  1.2];
zv = [1 1];
hold on
plot3(xv,yv,zv,'-.b');
grid on
%%%mcc setup
ai = analoginput('mcc');
addchannel(ai,0:2);
fs=64000; %升采样倍数
set(ai,'samplerate',fs);
T=0.2;      %s: 搜集数据T秒，计算一次位置
buffer = fs * T; %单次采样时间，对应采样点数
set(ai,'SamplesPerTrigger',buffer);
thta0 = 360;

% ii=0;
while(Gtrigger ==1)
    
%     ii=ii+1;
%     theta1 = (ii/40-fix(ii/40/2))*pi;
    if (Ttrigger ==0)
        pause(1);
    else
        theta =  Location_Own(fs,T,ai,thta0);
        thta0 = theta;
        theta1 = theta/180*pi;
        x1=cos(theta1);
        y1=sin(theta1);
        z1 = 1;
        h = findobj('color','red');
        delete(h);
        hold on
        plot3(x1,y1,z1,'or');   
    end
end

varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Ttrigger
Ttrigger =1


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Ttrigger
Ttrigger =0


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Gtrigger
Gtrigger =0
delete(gcbf);
