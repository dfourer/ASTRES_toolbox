function varargout = segtool(varargin)
% SEGTOOL MATLAB code for segtool.fig
%
%      SEGTOOL, by itself, creates a new SEGTOOL or raises the existing
%      singleton*.
%
%      H = SEGTOOL returns the handle to a new SEGTOOL or the handle to
%      the existing singleton*.
%
%      SEGTOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGTOOL.M with the given input arguments.
%
%      SEGTOOL('Property','Value',...) creates a new SEGTOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segtool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segtool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help segtool

% Last Modified by GUIDE v2.5 22-Feb-2012 11:03:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segtool_OpeningFcn, ...
                   'gui_OutputFcn',  @segtool_OutputFcn, ...
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


% --- Executes just before segtool is made visible.
function segtool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segtool (see VARARGIN)

% Choose default command line output for segtool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes segtool wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% Initialization
setappdata(handles.figure1,'mywav','cmor8-1');

% --- Outputs from this function are returned to the command line.
function varargout = segtool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --------------------------------------------------------------------
function uipushtool1_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Opens a signal contained in a MAT-file
[FileName,PathName,FilterIndex] = uigetfile({'*.mat','Matlab files (*.mat)'},'Please select a file');
mf = importdata([PathName FileName]);
s = mf.s;
dt = mf.dt;
if ~isvector(s) || ~isnumeric(dt)
    errordlg('Incorrect Matlab file');
    return;
end
noct = log2(length(s))-1;
if mod(noct,1) ~= 0
    warndlg('Warning, the size of the signal is not a power of 2! truncation of the signal');
    s = s(1:2^(floor(noct)+1));
end
    
% Save data
setappdata(handles.figure1,'dt',dt);
setappdata(handles.figure1,'snob',s);
setappdata(handles.figure1,'s',s);

% Plots signal
N = length(s);
t = dt*(0:N-1);
cla(handles.axes3);
plot(t,s,'b','Parent',handles.axes3);
set(handles.popupmenu6,'Value',1);
set(handles.slider1,'Visible','off');



% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load a predefined signal
sval = get(handles.popupmenu1,'Value');
N = 1024;
dt = 1/1023;
t = dt*(0:N-1);
switch(sval)
    case 1 % One sinus
        s = cos(2*pi*57*t);
    case 2 % Two sinus
        s = cos(2*pi*57*t) + 0.65*cos(2*pi*80*t);
    case 3 % One chirp
        s = cos(2*pi*27*(t+1.5*t.^2));
    case 4 % Two chirps
        s = cos(2*pi*27*(t+1.5*t.^2))+cos(2*pi*42*(t+1.5*t.^2));
    case 5 % sinus+chirp
        s = cos(2*pi*35*t)+cos(2*pi*45*(t+0.8*t.^2));
    case 6 % triangular signal
        s = sawtooth(20*pi*t)+1*sin(44*pi*t);
end

% Save data
setappdata(handles.figure1,'dt',dt);
setappdata(handles.figure1,'s',s);
setappdata(handles.figure1,'snob',s);

% Plots signal
cla(handles.axes3);
plot(t,s,'b','Parent',handles.axes3);
set(handles.popupmenu6,'Value',1);


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Launch the computation of synchrosqueezing

% get signal
s = getappdata(handles.figure1,'s');
dt = getappdata(handles.figure1,'dt');

% get wavelet/Fourier parameters
mywav = getappdata(handles.figure1,'mywav');

% get SQ parameters
nvval = [16 32 64];
nv = nvval(get(handles.popupmenu2,'Value'));
gamma = str2num(get(handles.edit2,'String'));

% Compute SQ
[Wx Tx fs as] = sqt(s,dt,gamma,mywav,nv);
Wx = flipud(Wx);
setappdata(handles.figure1,'Wx',Wx);

% WARNING!! to disappear
%figure();plot(abs(Wx(:,650)));xlabel('scale a');ylabel('|W(a,b_0)|');

setappdata(handles.figure1,'gamma',gamma);
setappdata(handles.figure1,'nv',nv);
setappdata(handles.figure1,'Tx',Tx);
setappdata(handles.figure1,'fs',fs);
setappdata(handles.figure1,'as',as);
sqplot(Wx,dt,as,handles.axes1);
sqplot(Tx,dt,fs,handles.axes2);
set(handles.popupmenu3,'Value',2);
set(handles.popupmenu4,'Value',3);


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

upplot(handles.figure1,handles.axes1,get(hObject,'Value'));


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu4.
function popupmenu4_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4

upplot(handles.figure1,handles.axes2,get(hObject,'Value'));


% --- Executes during object creation, after setting all properties.
function popupmenu4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3. Ridge extraction
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

dt = getappdata(handles.figure1,'dt');
mywav = getappdata(handles.figure1,'mywav');
nv = getappdata(handles.figure1,'nv'); 
Wx = getappdata(handles.figure1,'Wx');
as = getappdata(handles.figure1,'as');
reg=get(handles.checkbox1,'Value');

switch(get(handles.popupmenu8,'Value'))
    case 1 % Brevdo Ridge extraction
        % Get parameters
        Tx = getappdata(handles.figure1,'Tx');
        gamma = getappdata(handles.figure1,'gamma');
        fs = getappdata(handles.figure1,'fs');
        dt = getappdata(handles.figure1,'dt');
        nr = str2num(get(handles.edit4,'String'));
        lambda = str2num(get(handles.edit5,'String'));
        clwin = floor(str2num(get(handles.edit6,'String')));
        if isempty(Tx) ||isempty(fs) || isempty(dt) || isempty(mywav') || isempty(nv)
            warndlg('Some data is missing, extraction aborted');
            return;
        end
    
        % Ridge extraction
        [Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin);
        setappdata(handles.figure1,'Cs',Cs);
        setappdata(handles.figure1,'clwin',clwin);

        % Displaying Tx and Tx[Cs] together on another figure
        figure();
        ha = subplot(1,2,1);
        sqplot(Tx,dt,fs,ha);
        ha = subplot(1,2,2);
        sqplotRE(Tx,dt,fs,Cs,clwin,ha);
        
        % Computating a posteriori the wavelet mask SG
        SG = compute_mask(Wx,fs,Cs,clwin,dt,gamma);
        
        % Synthesis
        if reg==0
            imf = synth_tx(Tx,dt,mywav,nv,Cs,clwin);
        else
            imf = synth_mult(Wx,dt,mywav,nv,SG,1);
        end
        
    case 2 % Modified Brevdo 
        % Get parameters
        Tx = getappdata(handles.figure1,'Tx');
        gamma = getappdata(handles.figure1,'gamma');
        fs = getappdata(handles.figure1,'fs');
        dt = getappdata(handles.figure1,'dt');
        nr = str2num(get(handles.edit4,'String'));
        lambda = str2num(get(handles.edit5,'String'));
        clwin = floor(str2num(get(handles.edit6,'String')));
        if isempty(Tx) ||isempty(fs) || isempty(dt) || isempty(mywav') || isempty(nv)
            warndlg('Some data is missing, extraction aborted');
            return;
        end
    
        % Ridge extraction
        [Cs, Es] = brevridge_mult(Tx, fs, dt, nr, lambda, clwin);
        setappdata(handles.figure1,'Cs',Cs);
        setappdata(handles.figure1,'clwin',clwin);

        % Displaying Tx and Tx[Cs] together on another figure
        figure();
        ha = subplot(1,2,1);
        sqplot(Tx,dt,fs,ha);
        ha = subplot(1,2,2);
        sqplotRE(Tx,dt,fs,Cs,clwin,ha);
        
        % Computating a posteriori the wavelet mask SG
        Delta = floor(get_Delta(mywav)/(2^(1/nv)-1)+1-eps);
        centpsi = get_centpsi(mywav);
        SG = compute_maskWx(Wx,Cs,Delta,centpsi,fs);
        reg = get(handles.checkbox1,'Value');
        imf = synth_mult(Wx,dt,mywav,nv,SG,reg,as);
        
	case 3 % Meignen ridges
        % Get parameters
        s = getappdata(handles.figure1,'s');
        nmod = get(handles.edit4,'String');
        rrec = get(handles.edit5,'String');
        gamm = get(handles.edit6,'String');
        if isempty(Wx) ||isempty(as) || isempty(dt) || isempty(mywav)
            warndlg('Some data is missing, extraction aborted');
            return;
        end
    
        [SG Nf imf] = ts_segment(s,Wx,as,dt,mywav,nv,nmod,rrec,gamm,reg);

        
end
        
% Saving data
setappdata(handles.figure1,'SG',SG);
setappdata(handles.figure1,'imf',imf);

% Displaying imfs
imfplot(imf,dt,handles.axes4);

% Displaying Wx and SG together
sqplot(Wx,dt,as,handles.axes1);
set(handles.popupmenu3,'Value',3);
sqplot(exp(SG)-1,dt,as,handles.axes2);
set(handles.popupmenu4,'Value',4);

% Updating popupmunus 6 and 7 for displaying IMFs separately
nr = size(imf,1);
mycell = {'signal','Wavelet section','All IMFs','IMF 1','IMF 2','IMF 3','IMF 4',...
    'IMF 5','IMF 6','IMF 7','IMF 8','IMF 9','IMF 10','IMF 11','IMF 12',...
    'IMF 13','IMF 14','IMF 15','IMF 16','IMF 17','IMF 18','IMF 19','IMF 20',...
	'IMF 21','IMF 22','IMF 23','IMF 24','IMF 25','IMF 26','IMF 27','IMF 28'};
set(handles.popupmenu6,'String',mycell(1:3+nr));
set(handles.popupmenu7,'String',mycell(1:3+nr));
set(handles.popupmenu7,'Value',3);
set(handles.slider2,'Visible','off');


% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5

% Morlet/Morse wavelete
switch(get(hObject,'Value'))
    case 1 % Complex Morlet
        set(handles.text6,'String','Fc');
        set(handles.text7,'String','Fb');
        set(handles.edit7,'String','1');
        set(handles.edit8,'String','8');
    case 2 % Generalized Morse
        set(handles.text6,'String','beta');
        set(handles.text7,'String','gamma');
        set(handles.edit7,'String','7');
        set(handles.edit8,'String','3');
	case 3 % Bump Wavelet
        set(handles.text6,'String','mu');
        set(handles.text7,'String','sigma');
        set(handles.edit7,'String','1');
        set(handles.edit8,'String','0.2');
end


% --- Executes during object creation, after setting all properties.
function popupmenu5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4 : bouton Display in Wavelet
% panel
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Morlet/Morse wavelete
switch(get(handles.popupmenu5,'Value'))
    case 1 % Complex Morlet
        Fc = str2num(get(handles.edit7,'String'));
        Fb = str2num(get(handles.edit8,'String'));
        wavplot('cmor',Fb,Fc);
    case 2 % Generalized Morse
        beta = str2num(get(handles.edit7,'String'));
        gamma = str2num(get(handles.edit8,'String'));
        wavplot('gmor',beta,gamma);
	case 3 % Bump
        mu = str2num(get(handles.edit7,'String'));
        sigma = str2num(get(handles.edit8,'String'));
        wavplot('bump',mu,sigma);
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Change the mother wavelet
% Morlet/Morse wavelete
switch(get(handles.popupmenu5,'Value'))
    case 1 % Complex Morlet
        Fc = get(handles.edit7,'String');
        Fb = get(handles.edit8,'String');
        mywav = ['cmor' Fb '-' Fc];
    case 2 % Generalized Morse
        beta = get(handles.edit7,'String');
        gamma = get(handles.edit8,'String');
        mywav = ['gmor' beta '-' gamma];
	case 3 % Bump
        mu = get(handles.edit7,'String');
        sigma = get(handles.edit8,'String');
        mywav = ['bump' mu '-' sigma];
end
setappdata(handles.figure1,'mywav',mywav);
set(handles.uipanel2,'Title',['Mother wavelet ' mywav]);


% --------------------------------------------------------------------
function uipushtool2_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uipushtool2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Choosing which graphic is to save
prompt = {'Please choose the graphic to save : 1 (signal), 2 (IMFs), 3 4 or 5 (image plots of the wavelet transform/the Synchrosqueezed transform/the SQ with the curves))'};
dlg_title = 'Choosing the graphic';
num_lines = 1;
def = {'4'};
answer = inputdlg(prompt,dlg_title,num_lines,def);
if isempty(answer)
    return;
end


% "Saving", ie opening a new figure from axes
h = figure();
hga = axes('Parent',h);
switch(str2num(answer{1}))
    case 1
        ha = copyobj(handles.axes3,h);
    case 2
        ha = copyobj(handles.axes4,h);
    case 3
        ha = copyobj(handles.axes1,h);
    case 4
        ha = copyobj(handles.axes2,h);
end
set(ha,'Position',get(hga,'Position'));
delete(hga);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Wx = getappdata(handles.figure1,'Wx');
dt = getappdata(handles.figure1,'dt');
mywav = getappdata(handles.figure1,'mywav');
as = getappdata(handles.figure1,'as');
reconstruction(Wx,mywav,dt,as);


% --- Executes on selection change in popupmenu6.
function popupmenu6_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu6 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu6


changepopupv(hObject,handles.figure1,handles.axes3,handles.slider1);


% --- Executes during object creation, after setting all properties.
function popupmenu6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu7.
function popupmenu7_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu7 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu7

changepopupv(hObject,handles.figure1,handles.axes4,handles.slider2);

% --- Executes during object creation, after setting all properties.
function popupmenu7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slidercall(hObject,handles.figure1,handles.axes3);

% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

slidercall(hObject,handles.figure1,handles.axes4);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement : add Gaussian Noise
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

snob = getappdata(handles.figure1,'snob');
dt = getappdata(handles.figure1,'dt');
s = snob+norm(snob)/sqrt(length(snob))*get(hObject,'Value')*randn(size(snob));
setappdata(handles.figure1,'s',s);

% Plots signal
N = length(s);
t = dt*(0:N-1);
plot(t,s,'b','Parent',handles.axes3);
set(handles.popupmenu6,'Value',1);
set(handles.slider1,'Visible','off');




% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on selection change in popupmenu8. Choose the extraction
% method
function popupmenu8_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu8 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu8

update_param(hObject,handles);


% --- Executes during object creation, after setting all properties.
function popupmenu8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
