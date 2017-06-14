function varargout = reconstruction(varargin)
% RECONSTRUCTION MATLAB code for reconstruction.fig
%      RECONSTRUCTION, by itself, creates a new RECONSTRUCTION or raises the existing
%      singleton*.
%
%      H = RECONSTRUCTION returns the handle to a new RECONSTRUCTION or the handle to
%      the existing singleton*.
%
%      RECONSTRUCTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RECONSTRUCTION.M with the given input arguments.
%
%      RECONSTRUCTION('Property','Value',...) creates a new RECONSTRUCTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before reconstruction_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to reconstruction_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help reconstruction

% Last Modified by GUIDE v2.5 26-Jan-2012 14:46:03

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @reconstruction_OpeningFcn, ...
                   'gui_OutputFcn',  @reconstruction_OutputFcn, ...
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


% --- Executes just before reconstruction is made visible.
function reconstruction_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to reconstruction (see VARARGIN)

% Choose default command line output for reconstruction
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% Saving data
Wx = varargin{1};
mywav = varargin{2};
dt = varargin{3};
as = varargin{4};

setappdata(handles.figure1,'Wx',Wx);
setappdata(handles.figure1,'dt',dt);
setappdata(handles.figure1,'mywav',mywav);
setappdata(handles.figure1,'as',as);

% Creating graphics
% axes1 
N = size(Wx,2);
b = floor(N/2)+1;
plot(as,abs(Wx(:,b)),'Parent',handles.axes1);
set(handles.edit1,'String',num2str((b-1)*dt));
% axes2
tmp = 0.05*max(abs(Wx(:)));
set(handles.edit2,'String',num2str(tmp));
np = nb_peaks(tmp,Wx);
plot(dt*(0:N-1),np,'Parent',handles.axes2);
set(handles.text4,'String',num2str(mean(np)));
% Range of values for gamma
set(handles.edit3,'String',[num2str(0.2*tmp) '*[1:10]']);
% axes4
sqplot(Wx,dt,as,handles.axes4);


% UIWAIT makes reconstruction wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = reconstruction_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1 : plot vertical section  of
% Wx
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% getting the value t
t = str2num(get(handles.edit1,'String'));
dt = getappdata(handles.figure1,'dt');
Wx = getappdata(handles.figure1,'Wx');
as = getappdata(handles.figure1,'as');
N = size(Wx,2);
b = floor(t/dt)+1;
if b<1 || b>N
    warndlg(['Aborted : t must be a value between 0 and ' num2str(N*dt)]);
end
plot(as,abs(Wx(:,b)),'Parent',handles.axes1);



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


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


% --- Executes on button press in pushbutton2 : display graph 2
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

gamma = str2num(get(handles.edit2,'String'));
Wx = getappdata(handles.figure1,'Wx');
N = size(Wx,2);
dt = getappdata(handles.figure1,'dt');
np = nb_peaks(gamma,Wx);
plot(dt*(0:N-1),np,'Parent',handles.axes2);
set(handles.text4,'String',num2str(mean(np)));



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3 : compute ghat/gbar
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Getting parameters
eval(['Gam=' get(handles.edit3,'String') ';']);
if length(Gam)>50
    warndlg('To manue values in Gamma : restricting to only 20 values for computing considerations');
    Gam = linspace(Gam(1),Gam(end),50);
    set(handles.edit3,'String',[num2str(Gam(1)) ':' num2str(Gam(2)-Gam(1)) ':' num2str(Gam(end))]);
end
Wx = getappdata(handles.figure1,'Wx');
[gamma Nf flag] = comp_gamma(Wx,Gam);
set(handles.text8,'String',num2str(Nf));
set(handles.text9,'String',num2str(gamma));
if flag==0 % ghat
    set(handles.text7,'String','ghat');
else
    set(handles.text7,'String','gbar');
end
    
    


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


prompt = {'Please choose the graphic to save : 1, 2, 3, 4, or 5'};
dlg_title = 'Choosing the graphic';
num_lines = 1;
def = {'1'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return;
end


% "Saving", ie opening a new figure from axes
h = figure();
%set(h,'units','normalized');
hga = axes('Parent',h);
switch(str2num(answer{1}))
    case 1
        ha = copyobj(handles.axes1,h);
    case 2
        ha = copyobj(handles.axes2,h);
    case 3
        ha = copyobj(handles.axes3,h);
    case 4
        ha = copyobj(handles.axes4,h);
	case 5
        ha = copyobj(handles.axes5,h);
end
set(ha,'units','normalized');
set(ha,'Position',get(hga,'Position'));
delete(hga);
