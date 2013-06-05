function varargout = Gain(varargin)
% GAIN M-file for Gain.fig
%      GAIN, by itself, creates a new GAIN or raises the existing
%      singleton*.
%
%      H = GAIN returns the handle to a new GAIN or the handle to
%      the existing singleton*.
%
%      GAIN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GAIN.M with the given input arguments.
%
%      GAIN('Property','Value',...) creates a new GAIN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Gain_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Gain_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Gain

% Last Modified by GUIDE v2.5 30-Sep-2009 14:14:58

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Gain_OpeningFcn, ...
                   'gui_OutputFcn',  @Gain_OutputFcn, ...
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


% --- Executes just before Gain is made visible.
function Gain_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Gain (see VARARGIN)

% Choose default command line output for Gain
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Gain wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Gain_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

C=2.998e8;                                                          % (m/s)
NL_c=str2double(get(handles.Nonlinear_Parameter,'string'));         % (1/km/W)
bata_3=str2double(get(handles.Bata_3,'string'));                    % (ps^3/km)
bata_4=str2double(get(handles.Bata_4,'string'));                    % (ps^4/km)
wavelength_0=str2double(get(handles.Z_Wavelength,'string'));        % (nm)
wavelength_P=str2double(get(handles.P_Wavelength,'string'));        % (nm)
up=str2double(get(handles.UP,'string'));
down=str2double(get(handles.DOWN,'string'));
wavelength_S=linspace(down,up,(up-down)*2+1);                       % (nm)
Pp=str2double(get(handles.P_Power,'string'));                       % (W)

bata_2=bata_3*(C/wavelength_P-C/wavelength_0)*1e-3;                 % (ps^2/km)
det_fs=(C./wavelength_S-C/wavelength_P).*1e-3;                      % (1/ps)
det_bata=bata_2.*(det_fs.^2)+bata_4.*(det_fs.^4);                   % (1/km)
for k=1:(up-down)*2+1
    if (NL_c*Pp)^2-((det_bata(k)+2*NL_c*Pp)*0.5)^2>=0
        gain(k)=sqrt((NL_c*Pp)^2-((det_bata(k)+2*NL_c*Pp)*0.5)^2);	% (1/km)
    else
        gain(k)=0;
    end
end

% Drawing
axes(handles.axes1);
cla;
plot(wavelength_S,gain);
ylabel('Gain parameter (/km)');
xlabel('Signal Wavelength (nm)');
grid on;



function Nonlinear_Parameter_Callback(hObject, eventdata, handles)
% hObject    handle to Nonlinear_Parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Nonlinear_Parameter as text
%        str2double(get(hObject,'String')) returns contents of Nonlinear_Parameter as a double


% --- Executes during object creation, after setting all properties.
function Nonlinear_Parameter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Nonlinear_Parameter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_Power_Callback(hObject, eventdata, handles)
% hObject    handle to P_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_Power as text
%        str2double(get(hObject,'String')) returns contents of P_Power as a double


% --- Executes during object creation, after setting all properties.
function P_Power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_Power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function UP_Callback(hObject, eventdata, handles)
% hObject    handle to UP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of UP as text
%        str2double(get(hObject,'String')) returns contents of UP as a double


% --- Executes during object creation, after setting all properties.
function UP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to UP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function DOWN_Callback(hObject, eventdata, handles)
% hObject    handle to DOWN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of DOWN as text
%        str2double(get(hObject,'String')) returns contents of DOWN as a double


% --- Executes during object creation, after setting all properties.
function DOWN_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DOWN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function P_Wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to P_Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of P_Wavelength as text
%        str2double(get(hObject,'String')) returns contents of P_Wavelength as a double


% --- Executes during object creation, after setting all properties.
function P_Wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to P_Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bata_3_Callback(hObject, eventdata, handles)
% hObject    handle to Bata_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bata_3 as text
%        str2double(get(hObject,'String')) returns contents of Bata_3 as a double


% --- Executes during object creation, after setting all properties.
function Bata_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bata_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Bata_4_Callback(hObject, eventdata, handles)
% hObject    handle to Bata_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Bata_4 as text
%        str2double(get(hObject,'String')) returns contents of Bata_4 as a double


% --- Executes during object creation, after setting all properties.
function Bata_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Bata_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Z_Wavelength_Callback(hObject, eventdata, handles)
% hObject    handle to Z_Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Z_Wavelength as text
%        str2double(get(hObject,'String')) returns contents of Z_Wavelength as a double


% --- Executes during object creation, after setting all properties.
function Z_Wavelength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Z_Wavelength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


