function varargout = etGUI(varargin)

% etGUI will help you analyse the eye movements
% The core of the analysis is in showCurrentFrame()
% There is currently no documentation, but tooltips should be sufficient 
% to get you started using the GUI.
% 
% October 2014 - written by Michael Krumin


% ETGUI MATLAB code for etGUI.fig
%      ETGUI, by itself, creates a new ETGUI or raises the existing
%      singleton*.
%
%      H = ETGUI returns the handle to a new ETGUI or the handle to
%      the existing singleton*.
%
%      ETGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ETGUI.M with the given input arguments.
%
%      ETGUI('Property','Value',...) creates a new ETGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before etGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to etGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help etGUI

% Last Modified by GUIDE v2.5 26-Oct-2015 16:35:59

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @etGUI_OpeningFcn, ...
    'gui_OutputFcn',  @etGUI_OutputFcn, ...
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


% --- Executes just before etGUI is made visible.
function etGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to etGUI (see VARARGIN)

% Choose default command line output for etGUI
handles.output = hObject;

% default filepath
handles.filepath = '\\zserver.cortexlab.net\Data\EyeCamera\';
handles.filename = '';
handles.results = struct('x', [], 'y', [], 'aAxis', [], 'bAxis', [], 'abAxis', [], 'area', [], 'goodFit', [], 'blink', [],...
    'saturation', [], 'threshold', [], 'roi', []);
handles.doPlotting = true;
handles.state.roi = [];
handles.state.blinkROI = [];
handles.state.blinkThreshold = [];
handles.plotFigure = [];
handles.blinkPlot = [];

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes etGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = etGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in fastrewindTogglebutton.
function fastrewindTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to fastrewindTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

blinkplotState = get(handles.blinkplotPushbutton, 'Enable');
disableAll(handles.output);
set(hObject, 'Enable', 'on');
set(handles.frameSlider, 'Enable', 'inactive');
set(handles.iFrameText, 'Enable', 'on');
set(handles.frameSliderTitle, 'Enable', 'on');
set(handles.plotPushbutton, 'Enable', 'on');
set(handles.blinkplotPushbutton, 'Enable', blinkplotState);

while get(hObject, 'Value') & handles.iFrame>1
    handles.iFrame = handles.iFrame-1;
    guidata(hObject, handles);
    fastShowCurrentFrame(hObject, handles);
end

set(hObject, 'Value', 0);
enableAll(handles.output);
handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);


% --- Executes on button press in lastframePushbutton.
function lastframePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to lastframePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.iFrame>1)
    handles.iFrame = handles.iFrame-1;
    handles = showCurrentFrame(hObject, handles);
end

guidata(hObject, handles);

% --- Executes on button press in playTogglebutton.
function playTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to playTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

blinkplotState = get(handles.blinkplotPushbutton, 'Enable');
disableAll(handles.output);
set(hObject, 'Enable', 'on');
set(handles.plotPushbutton, 'Enable', 'on');
set(handles.blinkplotPushbutton, 'Enable', blinkplotState);
set(handles.saturationSlider, 'Enable', 'on');
set(handles.thresholdSlider, 'Enable', 'on');
set(handles.gaussEdit, 'Enable', 'on');
set(handles.diskEdit, 'Enable', 'on');
set(handles.originalCentreCheckbox, 'Enable', 'on');
set(handles.originalEdgeCheckbox, 'Enable', 'on');
set(handles.originalEllipseCheckbox, 'Enable', 'on');
set(handles.originalCropCheckbox, 'Enable', 'on');
set(handles.applyCheckbox, 'Enable', 'on');
set(handles.filteredCentreCheckbox, 'Enable', 'on');
set(handles.filteredEllipseCheckbox, 'Enable', 'on');
set(handles.saturationText, 'Enable', 'on');
set(handles.saturationTitle, 'Enable', 'on');
set(handles.thresholdText, 'Enable', 'on');
set(handles.thresholdTitle, 'Enable', 'on');
set(handles.gaussTitle, 'Enable', 'on');
set(handles.diskTitle, 'Enable', 'on');
set(handles.frameSlider, 'Enable', 'inactive');
set(handles.iFrameText, 'Enable', 'on');
set(handles.frameSliderTitle, 'Enable', 'on');
set(handles.text1, 'Enable', 'on');
set(handles.text2, 'Enable', 'on');


while get(hObject, 'Value') & handles.iFrame<handles.vr.NumberOfFrames
    handles.iFrame = handles.iFrame+1;
    handles = showCurrentFrame(hObject, handles);
    guidata(hObject, handles);
end

set(hObject, 'Value', 0);
enableAll(handles.output);

% --- Executes on button press in nextframePushbutton.
function nextframePushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to nextframePushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (handles.iFrame<handles.vr.NumberOfFrames)
    handles.iFrame = handles.iFrame+1;
    handles = showCurrentFrame(hObject, handles);
end

guidata(hObject, handles);

% --- Executes on button press in fastforwardTogglebutton.
function fastforwardTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to fastforwardTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

blinkplotState = get(handles.blinkplotPushbutton, 'Enable');
disableAll(handles.output);
set(hObject, 'Enable', 'on');
set(handles.frameSlider, 'Enable', 'inactive');
set(handles.iFrameText, 'Enable', 'on');
set(handles.frameSliderTitle, 'Enable', 'on');
set(handles.plotPushbutton, 'Enable', 'on');
set(handles.blinkplotPushbutton, 'Enable', blinkplotState);


while get(hObject, 'Value') & handles.iFrame<handles.vr.NumberOfFrames
    handles.iFrame = handles.iFrame+1;
    guidata(hObject, handles);
    fastShowCurrentFrame(hObject, handles);
end

set(hObject, 'Value', 0);
enableAll(handles.output);

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);


% --- Executes on slider movement.
function frameSlider_Callback(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

val = round(get(hObject, 'Value'));
set(hObject, 'Value', val);
handles.iFrame = val;

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);
% plotPushbutton_Callback(hObject, eventdata, handles);

% --- Executes during object creation, after setting all properties.
function frameSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to frameSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in roiPushbutton.
function roiPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to roiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disableAll(handles.output);

position = [handles.state.roi(1), handles.state.roi(2), diff(handles.state.roi(1:2:3)), diff(handles.state.roi(2:2:4))];
h = imrect(handles.originalAxis, position);
title(handles.originalAxis, 'Update the ROI and double-click')
newPos = wait(h);
handles.state.roi = [ceil(newPos(1)), ceil(newPos(2)), floor(newPos(1)+newPos(3)), floor(newPos(2)+newPos(4))];
delete(h);

enableAll(handles.output);

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);

% --- Executes on button press in originalEdgeCheckbox.
function originalEdgeCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to originalEdgeCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of originalEdgeCheckbox

showCurrentFrame(hObject, handles);

% --- Executes on button press in originalCentreCheckbox.
function originalCentreCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to originalCentreCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of originalCentreCheckbox

showCurrentFrame(hObject, handles);

% --- Executes on button press in originalEllipseCheckbox.
function originalEllipseCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to originalEllipseCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of originalEllipseCheckbox

showCurrentFrame(hObject, handles);

% --- Executes on slider movement.
function saturationSlider_Callback(hObject, eventdata, handles)
% hObject    handle to saturationSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% rounding the value to the second digit
val = get(hObject, 'Value');
val = round(val*100)/100;
set(hObject, 'Value', val);
handles.state.saturation = val;
handles = showCurrentFrame(hObject, handles);
satMode = '';
if ~isempty(handles.results.saturation) && ~all(isnan(handles.results.saturation))
    satMode = sprintf('   [%.2f]', mode(handles.results.saturation( ...
        ~isnan(handles.results.saturation))));
end
set(handles.saturationText, 'String', [num2str(get(hObject, 'Value')) satMode]);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function saturationSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saturationSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.state.saturation = get(hObject, 'Value');
guidata(hObject, handles);

% --- Executes on slider movement.
function thresholdSlider_Callback(hObject, eventdata, handles)
% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% rounding the value to the second digit
val = get(hObject, 'Value');
val = round(val*100)/100;
set(hObject, 'Value', val);
handles.state.threshold = val;
handles = showCurrentFrame(hObject, handles);
threshMode = '';
if ~isempty(handles.results.threshold) && ~all(isnan(handles.results.threshold))
    threshMode = sprintf('   [%.2f]', mode(handles.results.threshold( ...
        ~isnan(handles.results.threshold))));
end
set(handles.thresholdText, 'String', [num2str(get(hObject, 'Value')) threshMode]);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function thresholdSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to thresholdSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

handles.state.threshold = get(hObject, 'Value');
guidata(hObject, handles);


function gaussEdit_Callback(hObject, eventdata, handles)
% hObject    handle to gaussEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gaussEdit as text
%        str2double(get(hObject,'String')) returns contents of gaussEdit as a double

str = get(hObject, 'String');
num = str2num(str);
if ~isempty(num) && num>0
    handles.state.gauss = num;
    handles = showCurrentFrame(hObject, handles);
else
    set(hObject, 'String', num2str(handles.state.gauss));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function gaussEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gaussEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '0.5');
handles.state.gauss = 0.5;
guidata(hObject, handles);


function diskEdit_Callback(hObject, eventdata, handles)
% hObject    handle to diskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of diskEdit as text
%        str2double(get(hObject,'String')) returns contents of diskEdit as a double

str = get(hObject, 'String');
num = str2num(str);
if ~isempty(num) & num>0
    handles.state.disk = num;
    handles = showCurrentFrame(hObject, handles);
else
    set(hObject, 'String', num2str(handles.state.disk));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function diskEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to diskEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', '3');
handles.state.disk = 3;
guidata(hObject, handles);

% --- Executes on button press in filteredEdgeCheckbox.
function filteredEdgeCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to filteredEdgeCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filteredEdgeCheckbox


% --- Executes on button press in filteredCentreCheckbox.
function filteredCentreCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to filteredCentreCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filteredCentreCheckbox

showCurrentFrame(hObject, handles);

% --- Executes on button press in filteredEllipseCheckbox.
function filteredEllipseCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to filteredEllipseCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of filteredEllipseCheckbox

showCurrentFrame(hObject, handles);

% --- Executes on button press in originalCropCheckbox.
function originalCropCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to originalCropCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of originalCropCheckbox

showCurrentFrame(hObject, handles);

% --------------------------------------------------------------------
function openfile_pushtool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to openfile_pushtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isempty(eventdata) || ~isfield(eventdata, 'file')
    [filename, filepath] = uigetfile('*.*', 'Open eye-tracking video file...', handles.filepath);
    if isequal(filename, 0)
        fprintf('No files were selected\n');
        return;
    end
else
    [filepath, filename, ext] = fileparts(eventdata.file);
    filename = [filename, ext];
end

vr = VideoReader(fullfile(filepath, filename));
handles.filepath = filepath;
handles.filename = filename;
set(handles.filenameText, 'String', fullfile(filepath, filename));
nFrames = vr.NumberOfFrames;
nn = fieldnames(handles.results);
for iName = 1:length(nn)
    handles.results = setfield(handles.results, nn{iName}, nan(nFrames, 1));
end
handles.results.roi = nan(nFrames, 4);
handles.results.equation = cell(nFrames, 1);
handles.results.xxContour = cell(nFrames, 1);
handles.results.yyContour = cell(nFrames, 1);
handles.results.xxEllipse = cell(nFrames, 1);
handles.results.yyEllipse = cell(nFrames, 1);

set(handles.frameSlider, 'Value', 1);
set(handles.frameSlider, 'Min', 1, 'Max', nFrames);
% set(handles.frameSlider, 'SliderStep', [ceil(nFrames/100), ceil(nFrames/10)]);
% set(handles.iFrameText, 'String', sprintf('%d/%d', 1, nFrames))
% frame = read(vr, 1);
handles.vr = vr;
handles.iFrame = 1;
if isempty(handles.state.roi)
    handles.state.roi = [1 1 vr.Width, vr.Height];
    handles.state.blinkROI = handles.state.roi;
elseif handles.state.roi(3)>vr.Width || handles.state.roi(4)>vr.Height
    warning('ROI exceeds frame boundaries, resetting ROIs to full frame. File: %s', handles.filename);
    handles.state.roi = [1 1 vr.Width, vr.Height];
    handles.state.blinkROI = handles.state.roi;
else
    % do nothing, keep the old ROIs
end

% checking if there is already processed data
[~, name, ~] = fileparts(filename);
% if event.overwrite exists and is true, then overwrite existing results
overwrite = ~isempty(eventdata) && isfield(eventdata, 'overwrite') && eventdata.overwrite;
if exist(fullfile(filepath, [name, '_processed.mat']), 'file') && ~overwrite
    data = load(fullfile(filepath, [name, '_processed.mat']));
    handles.results = data.results;
    handles.state = data.state;
    
    set(handles.saturationSlider, 'Value', handles.state.saturation);
    set(handles.saturationText, 'String', ...
        [num2str(handles.state.saturation) ...
        sprintf('   [%.2f]', mode(handles.results.saturation))]);
    
    set(handles.thresholdSlider, 'Value', handles.state.threshold);
    set(handles.thresholdText, 'String', ...
        [num2str(handles.state.threshold) ...
        sprintf('   [%.2f]', mode(handles.results.threshold))]);
    
    set(handles.gaussEdit, 'String', num2str(handles.state.gauss));
    set(handles.diskEdit, 'String', num2str(handles.state.disk));
    if isfield(handles.state, 'blinkThreshold') && ~isempty(handles.state.blinkThreshold)
        set(handles.blinkthresholdEdit, 'Enable', 'on', 'String', num2str(handles.state.blinkThreshold));
        set(handles.blinkplotPushbutton, 'Enable', 'on');
    end
end

if any(isnan(handles.results.blink))
    set(handles.blinkplotPushbutton, 'Enable', 'off');
end

% pre-allocating space in the handle.movie variable
% bytesPerFrame = ffInfo.bytes;
% [~,sys] = memory;
% nFrames2LoadMax = floor(sys.PhysicalMemory.Available/bytesPerFrame*memFraction);


guidata(hObject, handles);
showCurrentFrame(hObject, handles);


function handles = showCurrentFrame(hObject, handles)

try
    frame = read(handles.vr, handles.iFrame);
    frame = frame(:,:,1);
catch
    fprintf('Failed loading a video frame\n');
    return;
end

frameCropped = frame(handles.state.roi(2):handles.state.roi(4), handles.state.roi(1):handles.state.roi(3));
gaussStd = str2num(get(handles.gaussEdit, 'String'));
diskR = str2num(get(handles.diskEdit, 'String'));
sat = get(handles.saturationSlider, 'Value');
th = get(handles.thresholdSlider, 'Value');

parsIn.diskR = diskR;
parsIn.gaussStd = gaussStd;
parsIn.sat = sat;
parsIn.th = th;

[parsOut, frameThresh] = analyseSingleFrame(frameCropped, parsIn);

x0 = parsOut.x0;
y0 = parsOut.y0;
aAxis = parsOut.aAxis;
bAxis = parsOut.bAxis;
abAxis = parsOut.abAxis;
area = parsOut.area;
equation = parsOut.equation; % not used, because of the xy-shifts related to ROI
xx = parsOut.xx;
yy = parsOut.yy;
good = parsOut.good;
xxEll = parsOut.xxEll;
yyEll = parsOut.yyEll;

if handles.doPlotting
    axes(handles.filteredAxis);
    hold off;
    imagesc(frameThresh);
    colormap gray;
    axis equal tight
    hold on;
    contour(frameThresh, [th, th], 'r:');
    axis off
end

xShift = handles.state.roi(1)-1;
yShift = handles.state.roi(2)-1;

if handles.doPlotting
    
    if get(handles.filteredCentreCheckbox, 'Value')
        plot(handles.filteredAxis, x0, y0, 'ro');
    end
    
    if get(handles.filteredEllipseCheckbox, 'Value')
%         str=sprintf('(x-%d)^2/%d+(y-%d)^2/%d+(x-%d)*(y-%d)/%d-%d', ...
%             x0, aAxis^2, y0, bAxis^2, x0, y0, abAxis, 1);
%         hh = ezplot(handles.filteredAxis, str, [1 handles.state.roi(3)-handles.state.roi(1) 1 handles.state.roi(4)-handles.state.roi(2)]);
        hh = plot(xxEll, yyEll);
        title('');
        set(hh, 'LineWidth', 0.5, 'LineStyle', '-', 'Color', 'g')
    end
    
    axes(handles.originalAxis);
    hold off;
    imagesc(frame);
    colormap gray;
    axis equal tight
    hold on;
    
    if get(handles.originalEdgeCheckbox, 'Value')
        plot(handles.originalAxis, xx+xShift, yy+yShift, 'r.')
    end
    
    if get(handles.originalCentreCheckbox, 'Value')
        plot(handles.originalAxis, x0+xShift, y0+yShift, 'ro');
    end
    
    if get(handles.originalEllipseCheckbox, 'Value')
        hh = plot(xxEll+xShift, yyEll+yShift);
        title('');
        set(hh, 'LineWidth', 0.5, 'LineStyle', '-', 'Color', 'g')        %     axis off;
    end
    
    if get(handles.originalCropCheckbox, 'Value')
        xlim(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)]);
        ylim(handles.originalAxis, [handles.state.roi(2), handles.state.roi(4)]);
    else
        plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(1)], [handles.state.roi(2), handles.state.roi(4)], 'b:');
        plot(handles.originalAxis, [handles.state.roi(3), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(4)], 'b:');
        plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(2)], 'b:');
        plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(4), handles.state.roi(4)], 'b:');
        if ~isequal(handles.state.roi, handles.state.blinkROI)
            plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(1)], [handles.state.blinkROI(2), handles.state.blinkROI(4)], 'y:');
            plot(handles.originalAxis, [handles.state.blinkROI(3), handles.state.blinkROI(3)], [handles.state.blinkROI(2), handles.state.blinkROI(4)], 'y:');
            plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(3)], [handles.state.blinkROI(2), handles.state.blinkROI(2)], 'y:');
            plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(3)], [handles.state.blinkROI(4), handles.state.blinkROI(4)], 'y:');
        end
    end
    
    if ~good
        plot(xlim, ylim, 'm', fliplr(xlim), ylim, 'm');
    end
    
    blink = handles.results.blink(handles.iFrame);
    if (~isnan(blink) && blink)
        xCoords = xlim;
        yCoords = ylim;
        text(xCoords(2), yCoords(1), 'B', 'Color', [1 0 0], 'FontSize', 20, ...
            'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
    end
    axis off;
    drawnow;
    
end

% updating the frame slider position and text
set(handles.iFrameText, 'String', sprintf('%d/%d', handles.iFrame, handles.vr.NumberOfFrames))
set(handles.frameSlider, 'Value', handles.iFrame);

% updating the results structure
% blink = handles.results.blink(handles.iFrame);
if get(handles.applyCheckbox, 'Value')
    if any(imag([x0, y0, aAxis, bAxis, area])) || ~good
        handles.results.x(handles.iFrame) = NaN;
        handles.results.y(handles.iFrame) = NaN;
        handles.results.aAxis(handles.iFrame) = NaN;
        handles.results.bAxis(handles.iFrame) = NaN;
        handles.results.abAxis(handles.iFrame) = NaN;
        handles.results.area(handles.iFrame) = NaN;
        handles.results.equation{handles.iFrame} = '';
        handles.results.xxContour{handles.iFrame} = [];
        handles.results.yyContour{handles.iFrame} = [];
        handles.results.xxEllipse{handles.iFrame} = [];
        handles.results.yyEllipse{handles.iFrame} = [];
    else
        handles.results.x(handles.iFrame) = x0+xShift;
        handles.results.y(handles.iFrame) = y0+yShift;
        handles.results.aAxis(handles.iFrame) = aAxis;
        handles.results.bAxis(handles.iFrame) = bAxis;
        handles.results.abAxis(handles.iFrame) = abAxis;
        handles.results.area(handles.iFrame) = area;
        handles.results.equation{handles.iFrame} = ...
            sprintf('(x-%d)^2/%d+(y-%d)^2/%d+(x-%d)*(y-%d)/%d-%d', ...
            x0+xShift, aAxis^2, y0+yShift, bAxis^2, x0+xShift, y0+yShift, abAxis, 1);
        handles.results.xxContour{handles.iFrame} = xx+xShift;
        handles.results.yyContour{handles.iFrame} = yy+yShift;
        handles.results.xxEllipse{handles.iFrame} = xxEll+xShift;
        handles.results.yyEllipse{handles.iFrame} = yyEll+yShift;
    end
    handles.results.goodFit(handles.iFrame) = good;
    handles.results.saturation(handles.iFrame) = get(handles.saturationSlider, 'Value');
    handles.results.threshold(handles.iFrame) = get(handles.thresholdSlider, 'Value');
    handles.results.roi(handles.iFrame, :) = handles.state.roi;
end

function fastShowCurrentFrame(hObject, handles)

frame = read(handles.vr, handles.iFrame);

axes(handles.originalAxis);
hold off;
imagesc(frame(:,:,1));
colormap gray;
axis equal tight
hold on;

if get(handles.originalCropCheckbox, 'Value')
    xlim(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)]);
    ylim(handles.originalAxis, [handles.state.roi(2), handles.state.roi(4)]);
else
    plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(1)], [handles.state.roi(2), handles.state.roi(4)], 'b:');
    plot(handles.originalAxis, [handles.state.roi(3), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(4)], 'b:');
    plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(2)], 'b:');
    plot(handles.originalAxis, [handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(4), handles.state.roi(4)], 'b:');
end

blink = handles.results.blink(handles.iFrame);
if (~isnan(blink) && blink)
    xCoords = xlim;
    yCoords = ylim;
    text(xCoords(2), yCoords(1), 'B', 'Color', [1 0 0], 'FontSize', 20, ...
        'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
end

axis off;
drawnow;

% updating the frame slider position and text
set(handles.iFrameText, 'String', sprintf('%d/%d', handles.iFrame, handles.vr.NumberOfFrames))
set(handles.frameSlider, 'Value', handles.iFrame);


% --- Executes on button press in plotPushbutton.
function plotPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to plotPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nFrames = handles.vr.NumberOfFrames;

if isempty(handles.plotFigure)
    screenSize = get(0, 'ScreenSize');
    handles.plotFigure = figure('Position', ...
        round([10 0.42*screenSize(4) screenSize(3)-20 0.5*screenSize(4)]));
    guidata(hObject, handles);
else
    figure(handles.plotFigure)
end
set(handles.plotFigure, 'WindowButtonDownFcn', 'clickAction;');
set(handles.plotFigure, 'UserData', handles);

ax = zeros(1,3);
blinks = ~isnan(handles.results.blink) & (handles.results.blink==true);
frameAxis = 1:nFrames;
ax(1) = subplot(3, 4, 1:3);
plot(1:nFrames, handles.results.x)
hold on
plot(frameAxis(blinks), handles.results.x(blinks), '.m');%, 'MarkerSize', 1);
plot([handles.iFrame, handles.iFrame], ylim, 'k:');
axis tight
ylabel('X [px]')
where = ylim;
where = (where-mean(where))*0.95+mean(where);
what = find(handles.results.goodFit==false);
plot(what, where(1)*ones(size(what)), '.m');
what = ones(size(frameAxis));
plot(frameAxis(blinks), where(2)*what(blinks), '.r');
ylim(prctile(handles.results.x, [0.1 99.9]));
hold off
title('Control-click to go to a specific frame in the main GUI');

ax(2) = subplot(3, 4, 5:7);
plot(1:nFrames, handles.results.y)
hold on
plot(frameAxis(blinks), handles.results.y(blinks), '.m');%, 'MarkerSize', 1);
plot([handles.iFrame, handles.iFrame], ylim, 'k:');
axis tight
ylabel('Y [px]')
where = ylim;
where = (where-mean(where))*0.95+mean(where);
what = find(handles.results.goodFit==false);
plot(what, where(1)*ones(size(what)), '.m');
what = ones(size(frameAxis));
plot(frameAxis(blinks), where(2)*what(blinks), '.r');
ylim(prctile(handles.results.y, [0.1 99.9]));
hold off

ax(3) = subplot(3, 4, 9:11);
plot(1:nFrames, handles.results.area)
hold on;
plot(frameAxis(blinks), handles.results.area(blinks), '.m');
plot([handles.iFrame, handles.iFrame], ylim, 'k:');
where = ylim;
where = (where-mean(where))*0.95+mean(where);
what = find(handles.results.goodFit==false);
plot(what, where(1)*ones(size(what)), '.m');
what = ones(size(frameAxis));
plot(frameAxis(blinks), where(2)*what(blinks), '.r');
hold off;
xlabel('Frame #');
ylabel('Area [px^2]')
axis tight
ylim(prctile(handles.results.area, [0.1 99.9]));


linkaxes(ax, 'x');

subplot(3,4,12)
frame = read(handles.vr, handles.iFrame);
hold off;
imagesc(frame(:,:,1));
colormap gray;
axis equal tight
xlimmem = xlim;
ylimmem = ylim;
hold on;
plot(handles.results.x(~blinks), handles.results.y(~blinks), 'y.');
plot(handles.results.x(blinks), handles.results.y(blinks), 'm.', 'MarkerSize', 1);
xlabel('X [px]');
ylabel('Y [px[');
set(gca, 'YDir', 'reverse');
axis off % comment out if you want axis labels and ticks

if get(handles.originalCropCheckbox, 'Value')
    xlim([handles.state.roi(1), handles.state.roi(3)]);
    ylim([handles.state.roi(2), handles.state.roi(4)]);
else
    plot([handles.state.roi(1), handles.state.roi(1)], [handles.state.roi(2), handles.state.roi(4)], 'w:');
    plot([handles.state.roi(3), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(4)], 'w:');
    plot([handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(2), handles.state.roi(2)], 'w:');
    plot([handles.state.roi(1), handles.state.roi(3)], [handles.state.roi(4), handles.state.roi(4)], 'w:');
    xlim(xlimmem);
    ylim(ylimmem);
end


% --- Executes on button press in analyzeTogglebutton.
function analyzeTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of analyzeTogglebutton

blinkplotState = get(handles.blinkplotPushbutton, 'Enable');

disableAll(handles.output);
set(hObject, 'Enable', 'on');
set(handles.frameSlider, 'Enable', 'inactive');
set(handles.iFrameText, 'Enable', 'on');
set(handles.frameSliderTitle, 'Enable', 'on');
set(handles.plotPushbutton, 'Enable', 'on');
set(handles.blinkplotPushbutton, 'Enable', blinkplotState);


if get(hObject, 'Value')
    frames2analyze = find(isnan(handles.results.x));
    handles.doPlotting = false;
    applyValue = get(handles.applyCheckbox, 'Value');
    set(handles.applyCheckbox, 'Value', 1);
    for iFrame = 1:length(frames2analyze)
        handles.iFrame = frames2analyze(iFrame);
        handles = showCurrentFrame(hObject, handles);
        %         if ~mod(iFrame, 100)
        drawnow %update
        %         end
        guidata(hObject, handles);
        if ~get(hObject, 'Value')
            break;
        end
    end
    handles.doPlotting = true;
    set(handles.applyCheckbox, 'Value', applyValue);
end

set(hObject, 'Value', 0);
enableAll(handles.output);

satMode = sprintf('   [%.2f]', mode(handles.results.saturation( ...
    ~isnan(handles.results.saturation))));
set(handles.saturationText, 'String', ...
    [num2str(get(handles.saturationSlider, 'Value')) satMode]);
threshMode = sprintf('   [%.2f]', mode(handles.results.threshold( ...
    ~isnan(handles.results.threshold))));
set(handles.thresholdText, 'String', ...
    [num2str(get(handles.thresholdSlider, 'Value')) threshMode]);

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);

% --- Executes on button press in applyCheckbox.
function applyCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to applyCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of applyCheckbox

grayColor = [0.941, 0.941, 0.941];
orangeColor = [1, 0.694, 0.392];
orangeColorDark = [0.871, 0.49, 0];
if get(hObject, 'Value')
    set(handles.playTogglebutton, 'BackgroundColor', orangeColorDark);
    set(handles.nextframePushbutton, 'BackgroundColor', orangeColorDark);
    set(handles.lastframePushbutton, 'BackgroundColor', orangeColorDark);
else
    set(handles.playTogglebutton, 'BackgroundColor', grayColor);
    set(handles.nextframePushbutton, 'BackgroundColor', grayColor);
    set(handles.lastframePushbutton, 'BackgroundColor', grayColor);
end
showCurrentFrame(hObject, handles);


% --- Executes on button press in blinkTogglebutton.
function blinkTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to blinkTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of blinkTogglebutton

if get(hObject, 'Value')
    disableAll(handles.output);
    %     set(hObject, 'Enable', 'on');
    set(handles.frameSlider, 'Enable', 'inactive');
    set(handles.iFrameText, 'Enable', 'on');
    set(handles.frameSliderTitle, 'Enable', 'on');
    
    nFrames = handles.vr.NumberOfFrames;
    if nFrames>1000
        frames = zeros(handles.vr.Height, handles.vr.Width, 1000);
        try
            idx = sort(randperm(nFrames, 1000));
        catch
            % hack for older versions of Matlab
            tmp = randperm(nFrames);
            idx = sort(tmp(1:1000));
        end
        for iFrame = 1:1000
            frame = read(handles.vr, idx(iFrame));
            frames(:,:,iFrame) = frame(:,:,1);
        end
    else
        frames = read(handles.vr);
    end
    
    averageROIFrame = mean(frames([handles.state.blinkROI(2):handles.state.blinkROI(4)], [handles.state.blinkROI(1):handles.state.blinkROI(3)], :), 3);
    
    rho = nan(nFrames, 1);
    for iFrame = 1:nFrames
        frame = read(handles.vr, iFrame);
        frameCropped = frame([handles.state.blinkROI(2):handles.state.blinkROI(4)], [handles.state.blinkROI(1):handles.state.blinkROI(3)], 1);
        tmp = corrcoef(single(frameCropped(:)), averageROIFrame(:));
        rho(iFrame) = tmp(2);
        
        set(handles.iFrameText, 'String', sprintf('%d/%d', iFrame, handles.vr.NumberOfFrames))
        set(handles.frameSlider, 'Value', iFrame);
        drawnow% update
    end
    % conservative detection
    if isempty(handles.state.blinkThreshold)
        handles.state.blinkThreshold = mean(rho)-4.5*std(rho);
    end
    tmp = (rho<=handles.state.blinkThreshold);
    % and then dilation (a +-100ms window around detected blinks will be marked as blinks)
    order = ceil(0.1*handles.vr.Framerate);
    tmp = filtfilt(ones(order+1, 1), 1, double(tmp(:)));
    handles.results.blink = tmp>0;
    handles.results.blinkRho = rho;
    handles.state.blinkROIAverageFrame = averageROIFrame;
    
    set(handles.blinkthresholdEdit, 'String', num2str(handles.state.blinkThreshold));
    %     set(handles.blinkthresholdEdit, 'Enable', 'on');
    %     set(handles.blinkplotPushbutton, 'Enable', 'on');
    
end

set(hObject, 'Value', 0);

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);

enableAll(handles.output);


% --------------------------------------------------------------------
function clearResults_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to clearResults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nFrames = handles.vr.NumberOfFrames;
nn = fieldnames(handles.results);
for iName = 1:length(nn)
    if isequal(nn{iName}, 'blink') || isequal(nn{iName}, 'blinkRho')
        continue;
    end
    handles.results = setfield(handles.results, nn{iName}, nan(nFrames, 1));
end
handles.results.roi = nan(nFrames, 4);
handles.results.equation = cell(nFrames, 1);
handles.results.xxContour = cell(nFrames, 1);
handles.results.yyContour = cell(nFrames, 1);

guidata(hObject, handles);


% --------------------------------------------------------------------
function save_uipushtool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to save_uipushtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[~, name, ~] = fileparts(handles.filename);
name2save = fullfile(handles.filepath, [name, '_processed.mat']);

data2save = struct;
data2save.results = handles.results;
data2save.state = handles.state;
save(name2save, '-struct', 'data2save');

% --------------------------------------------------------------------
function batch_uipushtool_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to batch_uipushtool (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% passing the handle to the etGUI gui to the etBatch
% this way etBatch will be able 'click' etGUI buttons 'remotely'

etBatch(get(get(hObject, 'Parent'), 'Parent'));

% --- Executes on button press in blinkroiPushbutton.
function blinkroiPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to blinkroiPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

disableAll(handles.output);

position = [handles.state.blinkROI(1), handles.state.blinkROI(2), diff(handles.state.blinkROI(1:2:3)), diff(handles.state.blinkROI(2:2:4))];
h = imrect(handles.originalAxis, position);
title(handles.originalAxis, 'Update the ROI and double-click')
newPos = wait(h);
handles.state.blinkROI = [ceil(newPos(1)), ceil(newPos(2)), floor(newPos(1)+newPos(3)), floor(newPos(2)+newPos(4))];
delete(h);

handles = showCurrentFrame(hObject, handles);
guidata(hObject, handles);

enableAll(handles.output);


% --- Executes on button press in blinkplotPushbutton.
function blinkplotPushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to blinkplotPushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

nFrames = handles.vr.NumberOfFrames;

if isempty(handles.blinkPlot)
    screenSize = get(0, 'ScreenSize');
    handles.blinkPlot = figure('Position', ...
        round([10 0.04*screenSize(4) screenSize(3)-20 0.3*screenSize(4)]));
    guidata(hObject, handles);
else
    figure(handles.blinkPlot)
end
set(handles.blinkPlot, 'WindowButtonDownFcn', 'clickAction;');
set(handles.blinkPlot, 'UserData', handles);

frameAxis = 1:nFrames;
plot(frameAxis, handles.results.blinkRho);
hold on;
plot(xlim, [handles.state.blinkThreshold handles.state.blinkThreshold], 'r:');
plot([handles.iFrame, handles.iFrame], ylim, 'k:')
realBlinks = handles.results.blinkRho<=handles.state.blinkThreshold;
filteredBlinks = (handles.results.blink & ~realBlinks);
plot(frameAxis(realBlinks), handles.results.blinkRho(realBlinks), 'r.');
plot(frameAxis(filteredBlinks), handles.results.blinkRho(filteredBlinks), 'g.');
xlim([1 nFrames]);
xlabel('Frame #');
ylabel('\rho(iFrame, averageFrame)');
title('Blink detection correlation');
hold off;
title('Control-click to go to a specific frame in the main GUI');


function blinkthresholdEdit_Callback(hObject, eventdata, handles)
% hObject    handle to blinkthresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of blinkthresholdEdit as text
%        str2double(get(hObject,'String')) returns contents of blinkthresholdEdit as a double

str = get(hObject, 'String');
num = str2num(str);
if ~isempty(num) && num>=-1 && num<=1
    handles.state.blinkThreshold = num;
    handles = showCurrentFrame(hObject, handles);
    tmp = (handles.results.blinkRho<=handles.state.blinkThreshold);
    % and then dilation (a +-100ms window around detected blinks will be marked as blinks)
    order = ceil(0.1*handles.vr.Framerate);
    tmp = filtfilt(ones(order+1, 1), 1, double(tmp(:)));
    handles.results.blink = tmp>0;
    handles = showCurrentFrame(hObject, handles);
else
    set(hObject, 'String', num2str(handles.state.blinkThreshold));
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function blinkthresholdEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to blinkthresholdEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in replayTogglebutton.
function replayTogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to replayTogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of replayTogglebutton

blinkplotState = get(handles.blinkplotPushbutton, 'Enable');
disableAll(handles.output);
set(hObject, 'Enable', 'on');
set(handles.frameSlider, 'Enable', 'inactive');
set(handles.iFrameText, 'Enable', 'on');
set(handles.frameSliderTitle, 'Enable', 'on');
set(handles.plotPushbutton, 'Enable', 'on');
set(handles.blinkplotPushbutton, 'Enable', blinkplotState);


while get(hObject, 'Value') && handles.iFrame<handles.vr.NumberOfFrames
    handles.iFrame = handles.iFrame+1;
    guidata(hObject, handles);
    replayCurrentFrame(hObject, handles);
end

set(hObject, 'Value', 0);
enableAll(handles.output);

% handles = replayCurrentFrame(hObject, handles);
guidata(hObject, handles);

%=======================================================================
function disableAll(h)

ch = get(h, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'off');
    end
end

hh = guidata(h);
ch = get(hh.blinkUipanel, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'off');
    end
end

ch = get(hh.uitoolbar1, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'off');
    end
end

%=======================================================================
function enableAll(h)

ch = get(h, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'on');
    end
end

hh = guidata(h);
ch = get(hh.blinkUipanel, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'on');
    end
end

ch = get(hh.uitoolbar1, 'Children');
for i=1:length(ch)
    if isprop(ch(i), 'Enable')
        set(ch(i), 'Enable', 'on');
    end
end

% set(hh.replayTogglebutton, 'Enable', 'inactive');
if isempty(hh.state.blinkThreshold)
    set(hh.blinkthresholdEdit, 'Enable', 'off');
    set(hh.blinkplotPushbutton, 'Enable', 'off');
end

%========================================================================

function replayCurrentFrame(hObject, handles)

frame = read(handles.vr, handles.iFrame);
frame = frame(:,:,1);

iFrame = handles.iFrame;
res = handles.results;

axes(handles.originalAxis);
hold off;
imagesc(frame);
colormap gray;
axis equal tight
hold on;

if get(handles.originalEdgeCheckbox, 'Value')
    plot(handles.originalAxis, res.xxContour{iFrame}, res.yyContour{iFrame}, 'r.')
end

if get(handles.originalCentreCheckbox, 'Value')
    plot(handles.originalAxis, res.x(iFrame), res.y(iFrame), 'ro');
end

if get(handles.originalEllipseCheckbox, 'Value')
    try
    hh = plot(handles.originalAxis, res.xxEllipse{iFrame}, res.yyEllipse{iFrame});
    title('');
    set(hh, 'LineWidth', 0.5, 'LineStyle', '-', 'Color', 'g')
    catch
        % apparently, there is no ellipse data
    end
    %     axis off;
end

if get(handles.originalCropCheckbox, 'Value')
    xlim(handles.originalAxis, [res.roi(iFrame, 1), res.roi(iFrame, 3)]);
    ylim(handles.originalAxis, [res.roi(iFrame, 2), res.roi(iFrame, 4)]);
else
    plot(handles.originalAxis, [res.roi(iFrame, 1), res.roi(iFrame, 1)], [res.roi(iFrame, 2), res.roi(iFrame, 4)], 'w:');
    plot(handles.originalAxis, [res.roi(iFrame, 3), res.roi(iFrame, 3)], [res.roi(iFrame, 2), res.roi(iFrame, 4)], 'w:');
    plot(handles.originalAxis, [res.roi(iFrame, 1), res.roi(iFrame, 3)], [res.roi(iFrame, 2), res.roi(iFrame, 2)], 'w:');
    plot(handles.originalAxis, [res.roi(iFrame, 1), res.roi(iFrame, 3)], [res.roi(iFrame, 4), res.roi(iFrame, 4)], 'w:');
    if ~isequal(res.roi(iFrame, :), handles.state.blinkROI)
        plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(1)], [handles.state.blinkROI(2), handles.state.blinkROI(4)], 'y:');
        plot(handles.originalAxis, [handles.state.blinkROI(3), handles.state.blinkROI(3)], [handles.state.blinkROI(2), handles.state.blinkROI(4)], 'y:');
        plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(3)], [handles.state.blinkROI(2), handles.state.blinkROI(2)], 'y:');
        plot(handles.originalAxis, [handles.state.blinkROI(1), handles.state.blinkROI(3)], [handles.state.blinkROI(4), handles.state.blinkROI(4)], 'y:');
    end
end

if ~isnan(res.goodFit(iFrame)) && ~res.goodFit(iFrame)
    plot(xlim, ylim, 'm', fliplr(xlim), ylim, 'm');
end

blink = res.blink(iFrame);
if (~isnan(blink) && blink)
    xCoords = xlim;
    yCoords = ylim;
    text(xCoords(2), yCoords(1), 'B', 'Color', [1 0 0], 'FontSize', 20, ...
        'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top');
end
axis off;
drawnow;


% updating the frame slider position and text
set(handles.iFrameText, 'String', sprintf('%d/%d', handles.iFrame, handles.vr.NumberOfFrames))
set(handles.frameSlider, 'Value', handles.iFrame);



function GotoFrameEdit_Callback(hObject, eventdata, handles)
% hObject    handle to GotoFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of GotoFrameEdit as text
%        str2double(get(hObject,'String')) returns contents of GotoFrameEdit as a double

input = get(hObject, 'String');
value = str2double(input);
if ~isnan(value)
    % means that the input was numerical
    nFrames = handles.vr.NumberOfFrames;
    % limiting the frameNumber to be within the movie frames
    frameNumber = min(max(1, round(value)), nFrames);
    handles.iFrame = frameNumber;
    handles = showCurrentFrame(hObject, handles);
    guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function GotoFrameEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to GotoFrameEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


