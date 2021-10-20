function varargout = LDT_Matlab(varargin)
%
%------------------Institut für Hochenergiephysik Wien---------------------
%--------------Oesterreichische Akademie der Wissenschaften----------------
%---------------------International Linear Collider------------------------
%
% Name of program:    LiC_Detector_Toy.m
% ILC coordinator:    Winfried Mitaroff
% Supervisor:         Meinhard Regler
% Authors:            Manfred Valentan, based on a program by R. Fruehwirth,
%                     with contributions by M. Regler. A part of the helix 
%                     tracking code was originally written by P. Billoir 
%                     for DELPHI
%
%--------------------------------------------------------------------------
%==========================================================================
% This program has been developed by the ILC group at the Institute of High 
% Energy Physics of the Austrian Academy of Sciences in Vienna, Austria. 
% It is distributed for scientific usage only. No parts of the code may be 
% passed on to third parties without consent of the authors. If results 
% obtained by using this program are presented publicly, please give the 
% appropriate credits to the authors.
%==========================================================================
%LDT_MATLAB M-file for LDT_Matlab.fig
%      LDT_MATLAB, by itself, creates a new LDT_MATLAB or raises the existing
%      singleton*.
%
%      H = LDT_MATLAB returns the handle to a new LDT_MATLAB or the handle to
%      the existing singleton*.
%
%      LDT_MATLAB('Property','Value',...) creates a new LDT_MATLAB using the
%      given property value pairs. Unrecognized properties are passed via
%      varargin to LDT_Matlab_OpeningFcn.  This calling syntax produces a
%      warning when there is an existing singleton*.
%
%      LDT_MATLAB('CALLBACK') and LDT_MATLAB('CALLBACK',hObject,...) call the
%      local function named CALLBACK in LDT_MATLAB.M with the given input
%      arguments.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help LDT_Matlab

% Last Modified by GUIDE v2.5 21-Apr-2008 09:17:28

% Begin initialization code - DO NOT EDIT
interrupt=0;
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @LDT_Matlab_OpeningFcn, ...
    'gui_OutputFcn',  @LDT_Matlab_OutputFcn, ...
    'gui_LayoutFcn',  [ ], ...
    'gui_Callback',   [ ]);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% Add subroutine directories to the MatLab search path
path(path,'subroutines');
path(path,'subroutines/gui');
path(path,'subroutines/outputs');
path(path,'subroutines/simulation');
path(path,'subroutines/tests');

global unit octave
unit=1;                        % unit=1   : computation in [mm]
                                % unit=10  : computation in [cm]
                                % unit=1000: computation in [m]
octave=0;
return

% --- Executes just before LDT_Matlab is made visible.
function LDT_Matlab_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output_control args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   unrecognized PropertyName/PropertyValue pairs from the
%            command line (see VARARGIN)

% Choose default command line output_control for LDT_Matlab
handles.output = hObject;
parameterfile='simulation_parameters_Matlab.txt';
handles.parameterfile=parameterfile;
%handles.fnameresults='results.mat';
set(handles.fnameresults,'String','results.mat')
guidata(hObject, handles);


% Update handles structure
guidata(hObject, handles);
set(hObject,'position',[10 27 180 30])

% UIWAIT makes LDT_Matlab wait for user response (see UIRESUME)
% uiwait(handles.figure1);
disp(' ');
disp('             __    ___   ______             ');
disp('            (  )  (   \ (__  __)            ');
disp('             )(__  ) > )   )(    v2.0       ');
disp('            (____)(___/   (__)              ');
disp('         The Vienna Fast Simulation Tool    ');
disp(' ');
disp('Ready!');
disp(' ');



% --- Outputs from this function are returned to the command line.
function varargout = LDT_Matlab_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output_control args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output_control from handles structure
varargout{1} = handles.output;
return



% --- Executes on selection change in barrelsetup.
function barrelsetup_Callback(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns setup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from setup
%contents = get(hObject,'String');
action=get(hObject,'Value');
geomdir='geom/';
switch action
    case 1 %Select barrel geometry
        list=get(handles.selected_barrel_geometry,'String');
        [fname,dirpath]=uigetfile([geomdir,'*.bgeom']);
        list{1}=[geomdir,fname];
        if any(fname~=0)
            set(handles.selected_barrel_geometry,'String',list)
            guidata(hObject, handles);
        else
            %errordlg('No barrel geometry has been selected')
            return
        end
    case 2 %Modify geometry
        %fname=get(handles.selected_barrel_geometry,'String');
        %fname=fname{1};
        [fname,dirpath]=uigetfile([geomdir,'*.bgeom']);
        if any(fname~=0)
            handles.filename=fname;
            guidata(hObject, handles);
        else
            return
        end
        %edit(handles.filename);
        edit([geomdir,fname])
    case 3 %New geometry
        %path(path,'geometry files/');
        edit geom/empty.bgeom;
    case 4
        list=get(handles.selected_barrel_geometry,'String');
        [fname,dirpath]=uigetfile([geomdir,'*.bgeom']);
        list{end+1}=[geomdir,fname];
        if any(fname~=0)
            %set(handles.selected_barrel_geometry,'String',Selected barrel geometry: ',barrelname})
            set(handles.selected_barrel_geometry,'String',list)
            guidata(hObject, handles);
        else
            %errordlg('No barrel geometry has been selected')
            return
        end
    case 5
        list=get(handles.selected_barrel_geometry,'String');
        if length(list)>1
            temp=list;
            clear list
            list=temp(1:end-1);
            set(handles.selected_barrel_geometry,'String',list)
            guidata(hObject, handles);
        end
end
return

% --- Executes on selection change in forwardsetup.
function forwardsetup_Callback(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns setup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from setup
%contents = get(hObject,'String');
action=get(hObject,'Value');
geomdir='geom/';
switch action
    case 1 %Select forward geometry
        list=get(handles.selected_forward_geometry,'String');
        [fname,dirpath]=uigetfile([geomdir,'*.fgeom']);
        list{1}=[geomdir,fname];
        if any(fname~=0)
            set(handles.selected_forward_geometry,'String',list)
            guidata(hObject, handles);
        else
            %errordlg('No barrel geometry has been selected')
            return
        end
    case 2 %Modify geometry
        %fname=get(handles.selected_barrel_geometry,'String');
        %fname=fname{1};
        [fname,dirpath]=uigetfile([geomdir,'*.fgeom']);
        if any(fname~=0)
            handles.filename=fname;
            guidata(hObject, handles);
        else
            return
        end
        %edit(handles.filename);
        edit([geomdir,fname])
    case 3 %New geometry
        %path(path,'geometry files/');
        edit geom/empty.fgeom;
    case 4 % add forward geometry to list
        list=get(handles.selected_forward_geometry,'String');
        [fname,dirpath]=uigetfile([geomdir,'*.fgeom']);
        list{end+1}=[geomdir,fname];
        if any(fname~=0)
            set(handles.selected_forward_geometry,'String',list)
            guidata(hObject, handles);
        else
            %errordlg('No barrel geometry has been selected')
            return
        end
    case 5 % remove forward geometry from list
        list=get(handles.selected_forward_geometry,'String');
        if length(list)>1
            temp=list;
            clear list
            list=temp(1:end-1);
            set(handles.selected_forward_geometry,'String',list)
            guidata(hObject, handles);
        end
end
return

function display2D_Callback(hObject, eventdata, handles)
%Display geometry
global Flags GeomVersion
global unit Bz SPR
global whandle mhandle fidlog
convf=0.299792458*1e-3*unit;    % conversion factor [Gev/c T^(-1) m^(-1)]
str=get(handles.selected_barrel_geometry,'String');
nameb=str;
if strcmp(nameb,'None')
    errordlg('No barrel geometry has been selected')
    return
end
str=get(handles.selected_forward_geometry,'String');
namef=str;
if strcmp(namef,'None')
    errordlg('No forward geometry has been selected')
    return
end
namep=handles.parameterfile;
mhandle=handles.messages;
whandle=handles.warnings;
%display=1;
%[SPR,N]=LDT_ReadParameters(namep);
%[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=LDT_ReadGeometry(nameb{1},namef{1});
%merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2);
%[theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
% angles to seperate forward,
% intermediate, barrel and rear
%varargout{1}=geometry3D(nameb{1},namef{1},VTX,SIT,TPC,SET,FM1,...
%    FM2,RM1,RM2,SPR,theta1,theta2,theta3,theta4,display);
LDT_display2D(nameb,namef);
% calls function 'geometry', if desired,
% which draws a sketch of the arrangement
%fclose(fidlog);
%h=LDT_main('display2D',nameb,namef,namep,messagehandle,warninghandle);
%handles.display=h;
%guidata(hObject, handles);

return

function display3D_Callback(hObject, eventdata, handles)
%Display geometry
global Flags GeomVersion
global unit Bz SPR
global whandle mhandle
convf=0.299792458*1e-3*unit;    % conversion factor [Gev/c T^(-1) m^(-1)]
str=get(handles.selected_barrel_geometry,'String');
nameb=str;
if strcmp(nameb,'None')
    errordlg('No barrel geometry has been selected')
    return
end
str=get(handles.selected_forward_geometry,'String');
namef=str;
if strcmp(namef,'None')
    errordlg('No forward geometry has been selected')
    return
end
namep=handles.parameterfile;
mhandle=handles.messages;
whandle=handles.warnings;
%display=2;
%[SPR,N]=LDT_ReadParameters(namep);
%[VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2]=LDT_ReadGeometry(nameb{1},namef{1});
%merging(VTX,SIT,TPC,SET,FM1,FM2,RM1,RM2);
%[theta1,theta2,theta3,theta4]=limitangles; % calculates the limiting theta
% angles to seperate forward,
% intermediate, barrel and rear
%varargout{1}=geometry3D(nameb{1},namef{1},VTX,SIT,TPC,SET,FM1,...
%    FM2,RM1,RM2,SPR,theta1,theta2,theta3,theta4,display);
LDT_display3D(nameb,namef);
% calls function 'geometry', if desired,
% which draws a sketch of the arrangement
%h=LDT_main('display3D',nameb,namef,namep,messagehandle,warninghandle);
%handles.display=h;
%guidata(hObject, handles);
return

% --- Executes on selection change in outputs.
function outputs_Callback(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns setup contents as cell array
%        contents{get(hObject,'Value')} returns selected item from setup
%contents = get(hObject,'String');
action=get(hObject,'Value');
%geomdir='geom/';
switch action
    case 1 %View log file
        edit LDT.log
    case 2 %Pulls histograms
        fnameresults=get(handles.fnameresults,'String');
        %histograms(action,fnameresults);
        LDT_pulls(fnameresults);
    case 3 %Residuals histograms
        fnameresults=get(handles.fnameresults,'String');
        %histograms(action,fnameresults);
        LDT_residuals(fnameresults);
    case 4 %Impacts histograms
        fnameresults=get(handles.fnameresults,'String');
        %impacts(fnameresults);
        LDT_impacts(fnameresults);
    case 5 % Curves
        fnameresults=get(handles.fnameresults,'String');
        LDT_curves(fnameresults);
    case 6 %Rave file
        [fname,dirpath]=uiputfile('*.txt');
        filename=[dirpath,fname];
        %rave(filename,handles.messages);
        fnameresults=get(handles.fnameresults,'String');
        LDT2rave(filename,fnameresults);
    case 7 %JAS3 file
        [fname,dirpath]=uiputfile('*.txt');
        filename=[dirpath,fname];
        fnameresults=get(handles.fnameresults,'String');
        LDT2jas3(filename,fnameresults);
end
return

% --- Executes during object creation, after setting all properties.
function barrelsetup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'fontsize',9)
return

% --- Executes during object creation, after setting all properties.
function forwardsetup_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'fontsize',9)
return

%
% --- Executes during object creation, after setting all properties.
function display2D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end
%set(hObject,'fontsize',9)
return
%}

% --- Executes during object creation, after setting all properties.
function display3D_CreateFcn(hObject, eventdata, handles)
% hObject    handle to setup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
%if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%    set(hObject,'BackgroundColor','white');
%end
%set(hObject,'fontsize',9)
return
%}

% --- Executes on selection change in setparameters.
function setparameters_Callback(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameters
%contents = get(hObject,'String');
% Edit simulation parameters
change_simulation_parameters(handles.parameterfile);

% --- Executes on selection change in run.
function run_Callback(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameters
%contents = get(hObject,'String');
% Run simulation
str=get(handles.selected_barrel_geometry,'String');
%nameb=str{1};
nameb=str;
str=get(handles.selected_forward_geometry,'String');
%namef=str{1};
namef=str;
namep=handles.parameterfile;
messagehandle=handles.messages;
warninghandle=handles.warnings;
if ~exist(nameb{1})
    errordlg('No barrel geometry has been selected')
    return
end
if ~exist(namef{1})
    errordlg('No forward geometry has been selected')
    return
end
%LDT_main('run simulation',nameb,namef,namep,messagehandle,warninghandle);
LDT_run(nameb,namef,messagehandle,warninghandle);
set(handles.fnameresults,'String','results.mat')
guidata(hObject, handles);
return

% --- Executes during object creation, after setting all properties.
function parameters_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'fontsize',9)
return

%{
% --- Executes on selection change in outputs.
function outputs_Callback(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns parameters contents as cell array
%        contents{get(hObject,'Value')} returns selected item from parameters
contents = get(hObject,'String');
action=get(hObject,'Value');
switch action
    case 1 % View log file
        edit LDT_Matlab.log
end
return
%}

% --- Executes during object creation, after setting all properties.
function outputs_CreateFcn(hObject, eventdata, handles)
% hObject    handle to parameters (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'fontsize',9)
return

% --- Executes on button press in stop.
function stop_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Save current geometry file names
str=get(handles.selected_barrel_geometry,'String');
nameb=str;
str=get(handles.selected_forward_geometry,'String');
namef=str;
save defaultgeometry nameb namef
close all
%close(handles.output)
return

% --- Executes on button press in saveresults.
function saveresults_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global fnameresults
save temp handles hObject
clear all
load temp handles hObject
[fname,dirpath]=uiputfile('*.mat');
str=get(handles.fnameresults,'String');
load(str);
if fname
    save(fname)
    %load temp
    str=get(handles.messages,'String');
    str{2}=['Simulation results written to'];
    str{3}=[dirpath,fname];
    str{4}=[];
    set(handles.messages,'String',str);drawnow;
    delete('temp.mat')
    set(handles.fnameresults,'String',fname)
    guidata(hObject, handles);
end
return

% --- Executes on button press in loadresults.
function loadresults_Callback(hObject, eventdata, handles)
% hObject    handle to stop (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global fnameresults
[fname,dirpath]=uigetfile('*.mat');
if fname
    set(handles.fnameresults,'String',fname)
    guidata(hObject, handles);
end
return

% --- Executes during object creation, after setting all properties.
function selected_barrel_geometry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_barrel_geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.selected_barrel_geometry=hObject;
set(hObject,'fontsize',7)
% Load default geometry
str=get(handles.selected_barrel_geometry,'String');
if exist('defaultgeometry.mat')
    load defaultgeometry nameb
    str=nameb;
else
    str{1}='None';
end
set(handles.selected_barrel_geometry,'String',str);
guidata(hObject, handles);
return

% --- Executes during object creation, after setting all properties.
function selected_forward_geometry_CreateFcn(hObject, eventdata, handles)
% hObject    handle to selected_barrel_geometry (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
handles.selected_forward_geometry=hObject;
set(hObject,'fontsize',7)
% Load default geometry
str=get(handles.selected_forward_geometry,'String');
if exist('defaultgeometry.mat')
    load defaultgeometry namef
    str=namef;
else
    str{1}='None';
end
set(handles.selected_forward_geometry,'String',str);
guidata(hObject, handles);
return

%{
% --- Executes on selection change in output.
function output_Callback(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns output contents as cell array
%        contents{get(hObject,'Value')} returns selected item from output
contents = get(hObject,'String');
action=contents{get(hObject,'Value')}
return
%}

% --- Executes during object creation, after setting all properties.
function output_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
set(hObject,'fontsize',9)
return

% --- Executes during object creation, after setting all properties.
function warnings_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.warnings=hObject;
guidata(hObject, handles);
set(hObject,'fontweight','bold')
set(hObject,'foregroundcolor',[0 0.5 0])
set(hObject,'fontsize',8)

% --- Executes during object creation, after setting all properties.
function messages_CreateFcn(hObject, eventdata, handles)
% hObject    handle to output (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
handles.messages=hObject;
guidata(hObject, handles);
set(hObject,'fontsize',8)




% --- Executes on selection change in popupmenu5.
function popupmenu5_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu5 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu5


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






% --- Executes on button press in cancel.
function cancel_Callback(hObject, eventdata, handles)
global interrupt
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
interrupt=1;
return





% --- Executes on button press in load_results.
function load_results_Callback(hObject, eventdata, handles)
% hObject    handle to load_results (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


