% Begin initialization code - DO NOT EDIT
function varargout = pix2wavCal(varargin)
% PIX2WAVCAL MATLAB code for pix2wavCal.fig
%      PIX2WAVCAL, by itself, creates a new PIX2WAVCAL or raises the existing
%      singleton*.
%
%      H = PIX2WAVCAL returns the handle to a new PIX2WAVCAL or the handle to
%      the existing singleton*.
%
%      PIX2WAVCAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PIX2WAVCAL.M with the given input arguments.
%
%      PIX2WAVCAL('Property','Value',...) creates a new PIX2WAVCAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before pix2wavCal_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to pix2wavCal_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help pix2wavCal

% Last Modified by GUIDE v2.5 15-Sep-2016 18:06:54

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @pix2wavCal_OpeningFcn, ...
    'gui_OutputFcn',  @pix2wavCal_OutputFcn, ...
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
end
% End initialization code - DO NOT EDIT
function varargout = pix2wavCal_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
try
    varargout{1} = handles.output;
catch

end
% The figure can be deleted now
%delete(gcf);
end

%% OPENING
function pix2wavCal_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to pix2wavCal (see VARARGIN)

% Choose default command line output for samplePicker3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

passedInHandles = varargin{1};
handles.wn = passedInHandles.wn;
handles.rootPath = passedInHandles.rootPath;
handles.rootPath = passedInHandles.rootPath;
handles.calSPEsPath = passedInHandles.calSPEsPath;
handles.savedCalsPath = passedInHandles.savedCalsPath;
calSPEsPath = passedInHandles.calSPEsPath; % for local func use below

pxpx = [0; 0; 143; 166; 264; 507; 624; 641];
wnpx = [620.9; 795.8; 1001.4; 1031.8; 1155.3; 1450.5; 1583.1; 1602.3];
handles.wnpx = wnpx;
wnpx = cellstr(num2str(wnpx));
pxpx = cellstr(num2str(pxpx));
% Update handles structure
guidata(hObject, handles);

% set up the button titles
for x=1:length(wnpx)
    field = ['checkbox',num2str(x)];
    set(handles.(field),'String',wnpx(x))
    
    % this rest is for debugging, so I don't have to manually set all the points
    % each time I want to test the functionality of this GUI
    field2 = ['wnedit',num2str(x)];
    set(handles.(field2),'String',pxpx(x))
    
end

cd(calSPEsPath);

calSPEsNames = dir('*.SPE'); % get the file info for every .SPE file
calSPEsNames = {calSPEsNames(~[calSPEsNames.isdir]).name}; % and extract the string names
rawSPEs = calSPEsNames'; % rotate the matrix bc I think better in columns for some reason
set(handles.calibListTable,'String',rawSPEs)

% default plot the first entry in the table
plotSpectra(hObject, handles, 1)

% plot the standard jpg for comparison
axes(handles.axes2)
cla
set(gca,'Tag','axes2')
I = imread('polystyrene','png');
imshow(I)

% UIWAIT makes pix2wavCal wait for user response (see UIRESUME)
uiwait(handles.figure1);
end

function plotSpectra(hObject, handles, n)
% 1. clear any previously plotted background-corrected data
axes(handles.PSbeadAxes)
cla
set(handles.PSbeadAxes,'Tag','PSbeadAxes')

% 2. grab the name of the selected bead
if n == 1
    contents = cellstr(get(handles.calibListTable,'String'));
    selectedCalSPE = contents{2};
    try
    set(handles.calibListTable,'Value',2)
    catch
    end
elseif n == 2
    contents = cellstr(get(hObject,'String'));
    selectedCalSPE = contents{get(hObject,'Value')};
end

data = double(readSPE(selectedCalSPE));  % read in the file with the readSPE function

% the following code handles the case where the SPE file contains several
% dimensions of spectra, i.e. non-accumulated scans. These would be
% stored in the 3rd dim of the readSPE output 'data1'
if size(data,3) > 1
    data = sum(squeeze(data),2)'; % manually 'accumulate' the scans
end

xDist = 1:1340;
plot(xDist,data,'Color','k','LineWidth',1);
hold on
set(gca,'xlim',[1 1340])
set(gca,'ButtonDownFcn',@(hObject,eventdata)pix2wavCal('PSbeadAxes_ButtonDownFcn',hObject,eventdata,guidata(hObject)))

xlimits = get(handles.PSbeadAxes,'xLim');
xMin = num2str(xlimits(1));
xMax = num2str(xlimits(2));

ylimits = get(handles.PSbeadAxes,'yLim');
yMin = num2str(ylimits(1));
yMax = num2str(ylimits(2));

set(handles.xMinEdit,'String',xMin);
set(handles.xMaxEdit,'String',xMax);
set(handles.yMinEdit,'String',yMin);
set(handles.yMaxEdit,'String',yMax);
end

%% plot selected spectrum from listbox
function calibListTable_Callback(hObject, eventdata, handles)
plotSpectra(hObject, handles, 2)

end

%% MANGE CHECKBOXES
function checkbox1_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox2_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox3_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox4_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox5_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox6_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox7_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox8_Callback(hObject, eventdata, handles)
checkbox_Callback_general(hObject, handles)
end
function checkbox_Callback_general(hObject, handles)

% get the number of the checkbox
str = strsplit(get(hObject,'Tag'),'box');
num = str{2};
field = ['set',num];

tf = get(hObject,'Value');  % return the state of the button
if tf == 1
    set(handles.(field),'Enable','on')
else
    set(handles.(field),'Enable','off')
end

end

%% SET TEXT BOXES
function wnedit1_Callback(hObject, eventdata, handles)
if checkValue(hObject) == 0
    set(hObject,'String',',<set>')
end
end
function wnedit2_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit3_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit4_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit5_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit6_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit7_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function wnedit8_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function tf = checkValue(hObject)
% this callback makes sure no value between n and l can be entered
n = 0; % change these limits if desired
l = 1340.1;
selVal = str2double(get(hObject,'String')); %returns contents of the object as a double

if selVal > l
    tf = 0;
elseif selVal < n
    tf = 0;
else
    tf = 1;
end
end

function set1_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set2_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set3_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set4_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set5_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set6_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set7_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function set8_Callback(hObject, eventdata, handles)
num = wnEdit_Callback_general(hObject,handles);
handles.selPt = num;
guidata(hObject,handles);
end
function num = wnEdit_Callback_general(hObject,handles)

% get the number of the checkbox
str = strsplit(get(hObject,'Tag'),'et');
num = str{2};

field = ['wnedit',num];

set(handles.(field),'String','');  % clear string on click

% get the number of the checkbox
str = strsplit(get(hObject,'Tag'),'et');
num = str{2};

end
function PSbeadAxes_ButtonDownFcn(hObject, eventdata, handles)

a = get(gca,'currentpoint');
val = num2str(a(1,1));
try
    field = ['wnedit',num2str(handles.selPt)];
    set(handles.(field),'String',val)
catch
    fprintf('Must select a button to <set>\n')
end

end

%% scale limits did change
function xMinEdit_Callback(hObject, eventdata, handles)
changeLims(handles)
end
function xMaxEdit_Callback(hObject, eventdata, handles)
changeLims(handles)
end
function yMinEdit_Callback(hObject, eventdata, handles)
changeLims(handles)
end
function yMaxEdit_Callback(hObject, eventdata, handles)
changeLims(handles)
end

function changeLims(handles)

xMin = str2double(get(handles.xMinEdit,'String'));
xMax = str2double(get(handles.xMaxEdit,'String'));
yMin = str2double(get(handles.yMinEdit,'String'));
yMax = str2double(get(handles.yMaxEdit,'String'));

if (xMax > xMin) && (yMax > yMin)
set(handles.PSbeadAxes,'xlim',[xMin xMax]);
set(handles.PSbeadAxes,'ylim',[yMin yMax]);
else
        % reset the limits
    xlimits = get(handles.PSbeadAxes,'xLim');
    xMin = num2str(xlimits(1));
    xMax = num2str(xlimits(2));
    
    ylimits = get(handles.PSbeadAxes,'yLim');
    yMin = num2str(ylimits(1));
    yMax = num2str(ylimits(2));
    
    set(handles.xMinEdit,'String',xMin);
    set(handles.xMaxEdit,'String',xMax);
    set(handles.yMinEdit,'String',yMin);
    set(handles.yMaxEdit,'String',yMax);
    % then drop an error message
    fprintf('the max limits must be greater than the min limits')

end
end

%% SAVE OR CANCEL
function calcWNbutton_Callback(hObject, eventdata, handles)

% initialize as a cell array so we can delete the empty rows with cellfun
% later
M = cell(8,2);
for x=1:8
    % check first that the checkbox next to the number is checked
    field = ['checkbox',num2str(x)]; % generate the appropriate handles id
    if get(handles.(field),'Value') == 1
        % if the box is enabled, grab the standard wn value
        M(x,1) = get(handles.(field),'String');
        % and the user entered value
        field2 = ['wnedit',num2str(x)];
        M(x,2) = cellstr(get(handles.(field2),'String'));
    end
end

% delete all the empty rows (unchecked boxes)
M(all(cellfun(@isempty,M),2), : ) = [];

% then covert to a matrix to perform issorted
M = cellfun(@str2double,M);

% check to make sure the values are sorted in ascending order before
% making the pix2wavCal
if issorted(M(:,2),'rows')
    set(handles.saveCalButton,'Enable','on')
    wn = polyval(polyfit(M(:,2),M(:,1),2),1:1340); % generate the wn by poly fittinghandles.wn = wn;
    handles.wn = wn;
    guidata(hObject, handles);
else
    set(handles.saveCalButton,'Enable','off')
    error('the set values are not in linear ascending order')
end
end

function saveCalButton_Callback(hObject, eventdata, handles)

S.wn = handles.wn;
dataArrayToPass(1) = {S};

x = get(handles.centerWnEdit,'String');
y = get(handles.gratingNameEdit,'String');

newName = ['centWN-',x,' gr-',y,'.mat'];

cd(handles.savedCalsPath);
save(newName,'dataArrayToPass')
guidata(hObject, handles);

% resume
uiresume(handles.figure1);
% and get back to the main window
delete(handles.figure1);

end

function cancelButton_Callback(hObject, eventdata, handles)

cd(handles.rootPath)

if isequal(get(handles.figure1, 'waitstatus'), 'waiting')
    % The GUI is still in UIWAIT, use UIRESUME
    uiresume(handles.figure1);
end

%  either way just close it
delete(handles.figure1);

end

%% HANGING CreateFcns
function calibListTable_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit1_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit2_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit5_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit6_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit7_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function wnedit8_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function xMinEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function xMaxEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function yMinEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function yMaxEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function gratingNameEdit_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function centerWnEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
