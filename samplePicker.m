%% BEGIN - DO NOT EDIT
function varargout = samplePicker(varargin)


%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @samplePicker_OpeningFcn, ...
    'gui_OutputFcn',  @samplePicker_OutputFcn, ...
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
end
function varargout = samplePicker_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% Get defau lt command line output from handles structure
varargout{1} = handles.output;
end
%% END - DO NOT EDIT

%% Menu Items
function FileMenu_Callback(~, ~, ~) %#ok<DEFNU>
end
function OpenMenuItem_Callback(~, ~, ~) %#ok<DEFNU>
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end
end
function PrintMenuItem_Callback(~, ~, handles) %#ok<DEFNU>
printdlg(handles.figure1)
end
function CloseMenuItem_Callback(~, ~, handles) %#ok<DEFNU>
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end
delete(handles.figure1)
end

%% --- Executes just before samplePicker is made visible.
function samplePicker_OpeningFcn(hObject, ~, handles, varargin)
% Choose default command line output for samplePicker
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

set(0, 'DefaultFigureRendererMode', 'manual')
set(0,'DefaultFigureRenderer','zbuffer')

set(handles.figure1,'toolbar','figure'); % add the figure toolbar to the window
set(handles.figure1,'Name','samplePicker');

%% set up Paths
rootPath = '/Users/randycarney/Developer/MATLAB/Noodle/';
handles.rootPath = rootPath;
handles.cdir = rootPath; % initialize the current directory handle (cdir) to the rootPath
cd(rootPath);
% first column is the handles.path name, second column is the actual path from the root
paths = {
    'rawDataSPEsPath','spectra';...
    'backgroundSPEsPath','spectra/backgrounds';...
    'calSPEsPath','spectra/calibration';...
    'subroutinesPath','subroutines';...
    'functionsPath','functions';...
    'exportFigfunctionsPath','functions/export_fig';...
    'savedStacksPath','saved/stacks';...
    'savedCalsPath','saved/pix2wavCals'
    'standardsPath','spectra/standards';
    };
    
for x=1:size(paths,1)
    field = paths{x,1}; % get var name
    handles.(field) = [rootPath paths{x,2}]; % set path name
    addpath(handles.(field))
end

guidata(hObject, handles); % and save

%% Initialize and save wavenumber and property matrices to GUIDATA struct properties
handles.colorMatrix = {'g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y'};
colors = {'Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen'};
for x=1:length(colors)
    colors2{x} = rgb(colors(x));
end
handles.colorMatrix2 = colors2;
handles.colors = colors;
handles.colorMatrixRGBtrip = {[0 1 0],[1 0 0],[0 0 1],[0 0 0],[1 0 1],[0 1 1],[1 1 0],[1 0 0],[0 1 0],[0 0 1],[0 0 0],[1 0 1],[0 1 1],[1 1 0]};

wn = pix2wavenum();
handles.wn = wn;

% get the user input crop (wn -> wnn) and splitFit wavenumbers
cropAbove = str2double(get(handles.cropAboveLabel,'String'));
cropBelow = str2double(get(handles.cropBelowLabel,'String'));
handles.cropAbove = cropAbove;
handles.cropBelow = cropBelow;
splitFitAt = str2double(get(handles.splitFitAtEditbox,'String'));
handles.splitFitAt = splitFitAt;
guidata(hObject, handles); % and save

% get the limits
wnn = checkNan(handles);
handles.wnn = wnn; % and save them
guidata(hObject, handles);

[rawDatasets, rawBackgrounds] = loadListboxes(handles);

% Pass the datset and background info to the appropriate GUIDATA handles
handles.rawDatasets = rawDatasets;
handles.rawBackgrounds = rawBackgrounds;
handles.rawDatasetsToPass = [];
guidata(hObject,handles); % and then save them

% set the default values for the backgroundListBox
set(handles.backgroundListbox,'Value', (1:2));
handles.selectedBackgrounds = get(handles.backgroundListbox, 'Value');
guidata(hObject,handles); % and then save them

backgroundListbox_Callback(hObject, x, handles);

% update status bar
set(handles.feedbackBarText,'String','Finished loading files.'); % status bar update

% UIWAIT makes samplePicker wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

function [rawDatasets, rawBackgrounds] = loadListboxes(handles)
%% Load the names of the base datasets

% Load all the PS bead cal files into the listbox
set(handles.feedbackBarText,'String','Loading all wn calibration files...'); % status bar update
cd(handles.savedCalsPath); % change directory to the rawData path
psCalNames = dir('*.mat'); % get the file info for every .mat file
psCalNames = {psCalNames(~[psCalNames.isdir]).name}; % and extract the string names
rawCalNames = psCalNames'; % rotate the matrix bc I think better in columns for some reason
if length(rawCalNames) < 1
    rawCalNames = '(no cals on file)';
end
set(handles.pixelWNcalListbox,'String',rawCalNames);
set(handles.pixelWNcalListbox,'Value',2);

% Load all the background files into the listbox
set(handles.feedbackBarText,'String','Loading all raw data files...'); % status bar update
cd(handles.backgroundSPEsPath); % change directory to the rawData path
backgroundSPEsNames = dir('*.SPE'); % get the file info for every .SPE file
backgroundSPEsNames = {backgroundSPEsNames(~[backgroundSPEsNames.isdir]).name}; % and extract the string names
rawBackgrounds = backgroundSPEsNames';
if length(rawBackgrounds) < 1
    rawBackgrounds = '(no backgrounds on file)';
end

% Load all raw data files
% by searching through the appropriate folders and parse filenames with format:
% 'datsetBaseFilename 1',...,'datsetBaseFilename n'
% where 'n' is an integer increasing by 1, the maximum being the total number of samples for a particular dataset
set(handles.feedbackBarText,'String','Loading all raw data files...'); % status bar update
cd(handles.rawDataSPEsPath); % change directory to the rawData path
rawDataSPEsNames = dir('*.SPE'); % get the file info for every .SPE file
rawDataSPEsNames = {rawDataSPEsNames(~[rawDataSPEsNames.isdir]).name}; % and extract the string names
rawDataSPEsNames = rawDataSPEsNames';

%switch back to the root directory
cd(handles.rootPath);

% Next, parse each name by splitting the 'datsetBaseFilename' from 'n'
datsetBaseFilenameNum = cell(2); % initialize a matrix to store the parsed info
for x=1:length(rawDataSPEsNames)
    splitName = strsplit(rawDataSPEsNames{x},{' ','.'}); % split by space character
    datsetBaseFilenameNum(x,1) = cellstr(splitName{1}); % take the name and put it in column 1
    datsetBaseFilenameNum(x,2) = cellstr(splitName{2});
end

% For that matrix, find all the unique dataset names and their occurences
baseNames = datsetBaseFilenameNum(:,1);
[uniqueNames, ~, s] = unique(baseNames);
occ = histc(s, 1:numel(uniqueNames));

% and load them into the rawDatasets matrix column 1 and 2
rawDatasets(:,1) = uniqueNames; % unique dataset names
rawDatasets(:,2) = cellstr(num2str(occ)); % and their number of samples

% update the GUI window listboxes with the datasets and samples
if length(rawDatasets) < 1
    rawDatasets = '(no datasets present in the directory)';
end
set(handles.datasetListbox,'String',rawDatasets(:,1));
set(handles.backgroundListbox,'String',rawBackgrounds);
set(handles.feedbackBarText,'String','Finished populating all listboxes.');

end

function wn = pix2wavenum()

% pixel to wavenumber calibration from polystyrene bead
% (calibration for gr - 600, cw - 1650)
wnpx = [1001.4 1031.8 1155.3 1450.5 1583.1 1602.3];
pxpx = [143 166 264 507 624 641];

% (calibration for gr - 600, cw - 1650)
% wnpx = [620.9 795.8 1001.4 1031.8 1155.3 1450.5 1583.1 1602.3];
% pxpx = [228 356 511 536 670.5 880 1000 1017];

wn = polyval(polyfit(pxpx,wnpx,2),1:1340); % generate the wn by poly fittinghandles.wn = wn;

end
function wnn = checkNan(handles)

cropAbove = handles.cropAbove;
cropBelow = handles.cropBelow;
wn = handles.wn;

% crop the data unless user-input is empty or the upper limit is smaller than the lower limit
if ~isnan(cropAbove) && ~isnan(cropBelow) && (cropBelow<cropAbove)
    wnn = wn(wn>=cropBelow & wn<=cropAbove);
else
    wnn = wn;
end

end

%% --- Dataset panel
function refreshDatasetButton_Callback(hObject, ~, handles) %#ok<DEFNU>

% [fileName,pathName] = uigetfile('*.xlsx','Select the base file');
% % get the first part of the file name up to the '.xlsx'
% addpath(pathName); % add to pathlist
% filename_split = strsplit(fileName,'.');
% set(handles.datasetListText,'String',filename_split{1});

[rawDatasets, rawBackgrounds] = loadListboxes(handles);

% Pass the datset and background info to the appropriate GUIDATA handles
handles.rawDatasets = rawDatasets;
handles.rawBackgrounds = rawBackgrounds;
guidata(hObject,handles); % and then save them

end
function changeDirButton_Callback(hObject, ~, handles) %#ok<DEFNU>

fn = uigetdir(handles.cdir,'Select spectra directory');
if fn ~= 0
    handles.cdir = fn;
    set(handles.datasetListText,'string',handles.cdir);
    handles.fileNames = dir(fullfile(handles.cdir,'*.SPE'));
    guidata(hObject, handles);
    
    if ~isempty(handles.fileNames)
        set(handles.datasetListbox,'enable','on');
        cd(handles.cdir)
        rawDatasets = datasetListBuilder(handles);
        handles.rawDatasets = rawDatasets;
        guidata(hObject, handles);
        
    else
        set(handles.datasetListbox,'enable','off');
    end
    cd(handles.rootPath)
end
end
function rawDatasets = datasetListBuilder(handles)

rawDataSPEsNames = dir('*.SPE'); % get the file info for every .SPE file
rawDataSPEsNames = {rawDataSPEsNames(~[rawDataSPEsNames.isdir]).name}; % and extract the string names
rawDataSPEsNames = rawDataSPEsNames';

% Next, parse each name by splitting the 'datsetBaseFilename' from 'n'
datsetBaseFilenameNum = cell(2); % initialize a matrix to store the parsed info
for x=1:length(rawDataSPEsNames);
    splitName = strsplit(rawDataSPEsNames{x},{' ','.'}); % split by space character
    datsetBaseFilenameNum(x,1) = cellstr(splitName{1}); % take the name and put it in column 1
    datsetBaseFilenameNum(x,2) = cellstr(splitName{2});
end

% For that matrix, find all the unique dataset names and their occurences
baseNames = datsetBaseFilenameNum(:,1);
[uniqueNames, ~, s] = unique(baseNames);
occ = histc(s, 1:numel(uniqueNames));

% and load them into the rawDatasets matrix column 1 and 2
rawDatasets(:,1) = uniqueNames; % unique dataset names
rawDatasets(:,2) = cellstr(num2str(occ)); % and their number of samples

% update the GUI window listboxes with the datasets and samples
if length(rawDatasets) < 1
    rawDatasets = '(no datasets present in the directory)';
end
set(handles.datasetListbox,'Value',1);
set(handles.datasetListbox,'String',rawDatasets(:,1));
end
function datasetListbox_Callback(hObject, ~, handles)

% only allows single selection at the moment
allData = plotRawSpectra(handles);
handles.data = allData;
guidata(hObject, handles);
cd(handles.rootPath)
end
function allData = plotRawSpectra(handles)
cd(handles.cdir)
% 1. clear any previously plotted background-corrected data
axes(handles.backgroundSpectraAxes)
cla
set(gca,'Tag','backgroundSpectraAxes')

% 2. get the name of the selected dataset
contents = cellstr(get(handles.datasetListbox,'String'));
selectedDataset = contents{get(handles.datasetListbox,'Value')};

% set the # samples text box
numSamples = str2double(handles.rawDatasets{get(handles.datasetListbox,'Value'),2});
set(handles.numSamplesLabel,'String',numSamples);

% 3. plot the raw spectra for the chosen dataset
axes(handles.rawSpectraAxes) % make active and clear the appropriate axes
cla
set(gca,'Tag','rawSpectraAxes')
hold on

% this loop looks for any SPE files in filepath (~/rawData) matching the string
% of the selected dataset and plots the associated spectra
allData = zeros(numSamples,1340);
textFilenamesAll = cell(numSamples,1);
for x=1:numSamples
    textFilename = [selectedDataset,' ',num2str(x),'.SPE']; % generate the expected filenames
    textFilenamesAll{x} = textFilename; % store each name in a matrix used to set the samples listbox after the loop
    
    data1 = double(readSPE(textFilename)); % read in the file with the readSPE function
    
    % the following code handles the case where the SPE file contains several
    % dimensions of spectra, i.e. non-accumulated scans. These would be
    % stored in the 3rd dim of the readSPE output 'data1'
    if size(data1,3) > 1
        data1 = sum(squeeze(data1),2)'; % manually 'accumulate' the scans
    end
    lgMult = str2double(get(handles.lagrangeMultEdit,'String'));
    data = whittaker_smoother(data1,lgMult);
    
    allData(x,:) = data;
    
    plot(handles.wn,data,'Color',handles.colorMatrix2{x});
    hold on
    
end

set(gca,'YLim',[-inf*0.8 inf*1.2]);
set(gca,'XLim',[-inf inf]);

set(handles.samplesListbox,'String',textFilenamesAll);

% Disable the 'Add to Stack Button' until backgrounds are applied
set(handles.AddToStackButton, 'Enable', 'off');
set(handles.applyBackgroundButton, 'Enable','on');

% update the feedback bar
newString = strcat('Plotted the raw spectra for dataset:',{' '},selectedDataset);
set(handles.feedbackBarText,'String',newString);

% reset the value of the samples list box
set(handles.samplesListbox,'Value',1);

% reset the selectedDatasetText label on the rawSpectraAxes
set(handles.selectedDatasetText,'String',selectedDataset)

end
function newCal_Callback(~, ~, handles) %#ok<DEFNU>

pix2wavCal(handles);

cd(handles.savedCalsPath); % change directory to the rawData path
psCalNames = dir('*.mat'); % get the file info for every .mat file
psCalNames = {psCalNames(~[psCalNames.isdir]).name}; % and extract the string names
rawCalNames = psCalNames'; % rotate the matrix bc I think better in columns for some reason
if length(rawCalNames) < 1
    rawCalNames = '(no cals on file)';
end
set(handles.pixelWNcalListbox,'String',rawCalNames);
set(handles.pixelWNcalListbox,'Value',2);

cd(handles.rootPath);

end
function pixelWNcalListbox_Callback(hObject, ~, handles) %#ok<DEFNU>

% get the chosen calibration filename
contents = cellstr(get(hObject,'String'));
chosenCalFilename = contents{get(hObject,'Value')};

% load the file
cd(handles.savedCalsPath);
S = load(chosenCalFilename);
cd(handles.rootPath); % cd back home

% store the data in a convenient array
S = S.dataArrayToPass{1};
handles.wn = S.wn;
guidata(hObject, handles);

wnn = checkNan(handles);
handles.wnn = wnn;
guidata(hObject, handles);

% update feedback bar
newString = ['Successfully changed calibration file to: ',chosenCalFilename];
set(handles.feedbackBarText,'String',newString); % status bar update

datasetListbox_Callback(hObject, 1, handles)

end

%% --- Background panel
function backgroundListbox_Callback(hObject, ~, handles)

selectedBackgrounds = get(handles.backgroundListbox,'Value');
sel = length(selectedBackgrounds);

handles.selectedBackgrounds = selectedBackgrounds;

% 1. get the user input crop values for wavenumber (wn -> wnn)
cropAbove = str2double(get(handles.cropAboveLabel,'String')); % update in handles according to box
cropBelow = str2double(get(handles.cropBelowLabel,'String'));
handles.cropAbove = cropAbove;
handles.cropBelow = cropBelow;
guidata(hObject, handles); % and save

% get the limits with the checkNan() function
wnn = checkNan(handles);
handles.wnn = wnn; % and save them
guidata(hObject, handles);

blanks = Func_bg_poly(handles);
handles.blanks = blanks;
guidata(hObject, handles); % save handles

trueBlanks = []; %#ok<NASGU> % initialize
if sel > 1
    trueBlanks = mean(blanks(1:sel,:));
else
    trueBlanks = blanks(1:sel,:);
end

axes(handles.selBackgroundAxes)
cla
set(gca,'Tag','selBackgroundAxes');
hold on
set(gca, 'XTickLabel',[]);
set(gca, 'YLim',[-inf*0.8 inf*1.2]);
plot(handles.wnn,trueBlanks,'Color','k')
vline(1450)

end
function autoFlOnlyCheckbox_Callback(hObject, ~, handles) %#ok<DEFNU>
tf = get(hObject,'Value');  % return the state of the button
if tf == 1
    set(handles.backgroundListbox,'Enable','off')
else
    set(handles.backgroundListbox,'Enable','on')
end
end
function splitFitCheckbox_Callback(hObject, ~, handles) %#ok<DEFNU>
tf = get(hObject,'Value');  % return the state of the button
if tf == 1
    set(handles.splitFitAtEditbox,'Enable','on')
    set(handles.splitFitAtEditbox,'String','1600')
else
    set(handles.splitFitAtEditbox,'Enable','off')
    set(handles.splitFitAtEditbox,'String','')
end
end
function splitFitAtEditbox_Callback(hObject, ~, handles) %#ok<DEFNU>

splitFitAt = str2double(get(hObject,'String'));

if isnan(splitFitAt) || (splitFitAt >= handles.cropAbove) || (splitFitAt <= handles.cropBelow)
    fprintf('Must selected a valid wavenumber at which to split the fit\n')
    set(hObject,'String',num2str(handles.splitFitAt));
    return
end

handles.splitFitAt = splitFitAt;
guidata(hObject, handles);end
function applyBackgroundButton_Callback(hObject, ~, handles) 
% housekeeping: initialize/clear the matrix where we will store the background corrected spectra
handles.data_mat = [];

% 1. get the user input crop values for wavenumber (wn -> wnn)
cropAbove = str2double(get(handles.cropAboveLabel,'String')); % update in handles according to box
cropBelow = str2double(get(handles.cropBelowLabel,'String'));
handles.cropAbove = cropAbove;
handles.cropBelow = cropBelow;
guidata(hObject, handles); % and save before getting the limits with the checkNan() function
wnn = checkNan(handles);
handles.wnn = wnn;
guidata(hObject, handles);

% 2. Use the selected backgrounds to compile a 'blanks' matrix
% make sure the user wants to background subtract anything at all
[blankMat, blankMat2] = Func_bg_poly(handles);
blanks = blankMat2;
bgRaw = mean(blankMat);

% 3. Prepare a few objects before graphing the BG-corrected graph
axes(handles.backgroundSpectraAxes) % make active the BG-corrected graph
cla % and clear it
set(gca,'Tag','backgroundSpectraAxes')
hold on

fields = {'data_mat','selectedDataset','sampleNum','sampleColors'};
for x=1:length(fields)
    try
        handles = rmfield(handles,fields{x});
    catch
    end
end

% 4. Iterate through each sample, correct with blanks matrix, and plot
[numOfSamples, ~] = size(handles.data); % get the current number of samples
for x=1:numOfSamples
    data = handles.data(x,:);
    % Perform the background correction
    % split the fit if the user input checkbox is checked
    if get(handles.splitFitCheckbox,'Value') ~= 1
        data = data(handles.wn>=cropBelow & handles.wn<=cropAbove); % and crop it according to the user-input
        coeffs = asls(blanks',data',.05,[]); % estimate the background coefficients using AsLS
        corrected = data-(coeffs'*blanks);  % construct the estimated background and subtract it from data
        if x == 1
            bgEst = (coeffs'*blanks);
        end
    else
        splitFitAt = str2double(get(handles.splitFitAtEditbox,'String'));
        handles.splitFitAt = splitFitAt; 
        guidata(hObject, handles);
        data1 = data(handles.wn>=cropBelow & handles.wn<=splitFitAt);
        data2 = data(handles.wn>=splitFitAt & handles.wn<=cropAbove);
        for z=1:size(blanks,1)
            blanksRow = blanks(z,:);
            blanks1(z,:) = blanksRow(handles.wnn<=splitFitAt);
            blanks2(z,:) = blanksRow(handles.wnn>=splitFitAt);
        end
        coeffs1 = asls(blanks1',data1',.05,[]); % estimate the background coefficients using AsLS
        coeffs2 = asls(blanks2',data2',.05,[]);
        corrected1 = data1-(coeffs1'*blanks1);  % construct the estimated background and subtract it from data
        corrected2 = data2-(coeffs2'*blanks2);
        corrected = horzcat(corrected1,corrected2);
        if x == 1
            bgEst1 = (coeffs1'*blanks1);
            bgEst2 = (coeffs2'*blanks2);
            bgEst = horzcat(bgEst1,bgEst2);
        end
    end
    
    %normalization
    
    normrange = find(handles.wnn>cropBelow & handles.wnn<cropAbove);
    normarea = trapz(handles.wnn(normrange),corrected(:,normrange));
    normspec = corrected./repmat(normarea,1,size(corrected,2));
   
    % we want to set the color of the graph according the sample number N
    currentSampleNamesInGraph = get(handles.samplesListbox,'String'); % return the sample names
    selectedSampleName = strsplit(currentSampleNamesInGraph{x},{' ','.'}); % get the name of the selected sample
    selectedSampleN = str2double(selectedSampleName{2});% grab the sample number N
   
    if get(handles.normalizeDataCheckbox,'Value') == 1
        % plot the line with color matching the selectedSampleN
        plot(handles.wnn, normspec,'Color', handles.colorMatrix2{selectedSampleN})
        hold on
        handles.data_mat(x,:) = normspec;
    else
        plot(handles.wnn, corrected,'Color', handles.colorMatrix2{selectedSampleN})
        hold on
        handles.data_mat(x,:) = corrected;
    end
    
    % concatenate data for passing to Stack
    
    handles.selectedDataset(x,:) = selectedSampleName;
    handles.sampleNum(x,:) = x;
    handles.sampleColors(x) = handles.colors(x);
end
hline(0)

% plot the background example for the first spectra in the rawSpectraAxes
axes(handles.rawSpectraAxes) % make active the BG-corrected graph
% if there is a line there already, delete it
h = findobj(gcf,'Color',rgb('Green'),'LineStyle','--','LineWidth',1);
g = findobj(gcf,'Color',rgb('Green'),'LineStyle','-.','LineWidth',0.75);
if ~isempty(h)
    delete(h)
end
if ~isempty(g)
    delete(g)
end
hold on
plot(handles.wnn, bgEst,'Color',rgb('Green'),'LineStyle','--','LineWidth',1)
if get(handles.plotBGAvgCheckbox,'Value') == 1
    plot(handles.wnn, bgRaw,'Color',rgb('Green'),'LineStyle','-.','LineWidth',0.75)
end

set(handles.backgroundSpectraAxes, 'xlim',[cropBelow cropAbove]);
set(handles.rawSpectraAxes, 'xlim',[cropBelow cropAbove]);
% 5. save the handles.data_mat with all of the BG-corrected data
guidata(hObject, handles);

% 6. Housekeeping
% active the 'Add Dataset To Stack' Button
set(handles.AddToStackButton, 'Enable', 'on');

set(handles.selectedDatasetText,'String',selectedSampleName{1})

% update the feedback bar
newString = strcat('Plotted the background corrected dataset: ',selectedSampleName{1});
set(handles.feedbackBarText,'String',newString);

end
function [blankMat, blankMat2] = Func_bg_poly(handles)
% inputs the names of all backgrounds on file and returns a matrix of NxM
% background spectra, where N is # selected bgs and M is the value of each
% point in the spectra
allBackgroundNames = handles.rawBackgrounds;
backgroundSelection = handles.selectedBackgrounds;
wn = handles.wn;
cropAbove = handles.cropAbove;
cropBelow = handles.cropBelow;

% initialize the holding matrices for the upcoming 'for' loop
selectedBgNames = cell(size(backgroundSelection));

% iterate over each user-highlighted background spectra
num_backgrounds = length(backgroundSelection);

if get(handles.autoFlOnlyCheckbox,'Value') == 1
    blankSpectra = zeros(1,1340);
    blankMat = blankSpectra(wn>=cropBelow & wn<=cropAbove);
    autoFl = str2double(get(handles.autoFledit,'String'));
    % this matrix is output to the return of the called function
    blankMat2 = mypoly(autoFl,length(blankMat));
    return
end

for N = 1:num_backgrounds
    % 1. fill a matrix with the background names
    % (currently unused at the moment)
    selectedBgNames(1,N) = cellstr(allBackgroundNames{backgroundSelection(1,N)}); % create a cell array of bg names
    % 2. read in the current background SPE file in blankPlaceholder
    blankSpectra = double(readSPE(allBackgroundNames{backgroundSelection(1,N)})); % read each bg SPE into a matrix
    % 3. manually 'accumulate' the scans if the user did not specificy
    % accumulations during aquisition
    if size(blankSpectra,3) > 1
        blankSpectra = sum(squeeze(blankSpectra),2)';
    end
    % blankSpectra should now be of dimension: [1 x 1340]
    % 4. It is cropped according to the user-input values (passed in
    % handles)
    croppedBlankSpectra = blankSpectra(wn>=cropBelow & wn<=cropAbove);
    % 5. And finally smoothed
    % grab the Lagrange parameter from the user input box
    bgLG = str2double(get(handles.lagrangeMultiEditBG,'String'));
    qz1 = whittaker_smoother(croppedBlankSpectra,bgLG);
    
    % and appended to a loopholder matrix that gets passed back as this
    % function's output
    blankMat(N,:) = qz1;    % append the smoothed, chopped blank to a matrix
    
end

% finish constructing the background by appending a n-th order
% polynomial function representing system autofluorescence (will be used to
% fit to the data with asls)
% polynomial order is specified by user input:
autoFl = str2double(get(handles.autoFledit,'String'));

% this matrix is output to the return of the called function
blankMat2 = [blankMat;mypoly(autoFl,length(blankMat))];

end

%% --- Samples panel
function samplesListbox_Callback(hObject, ~, handles) %#ok<DEFNU>
% 1. when a slection is made, we want to store the selection in handles
selectedIndex = get(hObject,'Value');
handles.selectedIndex = selectedIndex;
guidata(hObject, handles);

% 2. then find the corresponding lines in the graphs and bold them
% The selected index may not match the N of the sample if previous lines
% were already deleted, so we have to find the 'N' for the selected sample
% a. find the N of the selected sample
currentSampleNamesInGraph = get(handles.samplesListbox,'String'); % get the list of names
try
    selectedSampleName = strsplit(currentSampleNamesInGraph{selectedIndex},{' ','.'}); % get the name of the selected sample
catch
    return
end
selectedSampleN =str2double(selectedSampleName{2});

% b. get the color that this sample would have been graphed in
color = handles.colorMatrix2{selectedSampleN};

% c. access the figure and grab the appropriate colored line in graph
h = findobj(0,'Tag','rawSpectraAxes');
h2 = findobj(0,'Tag','backgroundSpectraAxes');
% first remove any bolded lines from previous selections
allLines = findall(h,'type','line'); % return all graphed lines in both graphs
set(allLines,'LineWidth',0.5); % set all of ther widths to not-bold (0.5)
allLines = findall(h2,'type','line'); % return all graphed lines in both graphs
set(allLines,'LineWidth',0.5); % set all of ther widths to not-bold (0.5)
% then find the lines we want to bold
allLinesBG = findall(h2,'type','line'); % return the lines in the background graph
allLinesRawData = findall(h,'type','line'); % and raw data graph
selectedBGLine = findall(allLinesBG, 'Color',color); % and grab the ones matching the color we are looking for
selectedRawLine = findall(allLinesRawData, 'Color',color);
% Bold the selected lines
set(selectedBGLine,'LineWidth',2);
set(selectedRawLine,'LineWidth',2);

% 3. Finally store the selected lines in the GUIDATA handles for deletion if
% necessary
handles.selectedBGLine = selectedBGLine;
handles.selectedRawLine = selectedRawLine;
guidata(hObject, handles);

end

%% Stack panel
function AddToStackButton_Callback(hObject, ~, handles) %#ok<DEFNU>

try
    numInStack = size(handles.stack,1); % return the current length of the stack
catch
    numInStack = 0;
end
samplesListBox = get(handles.samplesListbox,'String');
numSamplesToAdd = size(samplesListBox,1); % and the number of samples to add

% set the range for which to add the data (below any previous data)
from = numInStack+1;
to = numInStack + numSamplesToAdd;

% then save into the proper range of the GUIDATA handles
handles.stack(from:to,:) = handles.data_mat(1:numSamplesToAdd,:);
handles.stackIDs(from:to,:) = handles.selectedDataset(1:numSamplesToAdd,:);
guidata(hObject, handles); % and save

% finally, update the Stack list in the stackListbox
a = get(handles.stackListbox,'String');
try
    a = vertcat(a,handles.selectedDataset(1,1));
catch
end
toAdd = [handles.selectedDataset(1,1),cellstr(num2str(size(handles.selectedDataset,1))),handles.colors{size(handles.rawDatasetsToPass,1)+1}];
handles.rawDatasetsToPass = [handles.rawDatasetsToPass; toAdd];
guidata(hObject, handles); % and save

set(handles.stackListbox,'String',a);
set(handles.stackListbox,'Value',length(a));

try
    newString = strcat('Added',{' '},handles.selectedDataset{1,1},{' '},'to Stack');
catch
    newString = 'No samples in dataset to add';
end

% clear the graphs and listboxes, enable 'Add dataset' and 'Apply PCA'

axes(handles.backgroundSpectraAxes)
cla
set(gca,'Tag','backgroundSpectraAxes');
axes(handles.rawSpectraAxes)
cla
set(gca,'Tag','rawSpectraAxes');
set(handles.cropAboveLabel, 'Enable','off');
set(handles.cropBelowLabel, 'Enable','off');
set(handles.samplesListbox,'String','');
set(handles.applyBackgroundButton, 'Enable','off');
set(handles.pixelWNcalListbox, 'Enable','off');
set(handles.AddToStackButton, 'Enable','off');
set(handles.applyPCAbutton, 'Enable','on');
guidata(hObject, handles); % and save

set(handles.feedbackBarText,'String',newString);

% automatically load up the raw spectra and BG for the next dataset in the block
str = get(handles.datasetListbox,'String');
val = get(handles.datasetListbox,'Value');

if val == size(str,1)
    set(handles.datasetListbox,'Value',val)
else
    set(handles.datasetListbox,'Value',(val+1))
end
allData = plotRawSpectra(handles);
handles.data = allData;
guidata(hObject, handles);
cd(handles.rootPath)
applyBackgroundButton_Callback(hObject, 1, handles) % plot backgrrounds


end

%% DATA HANDLING (STACK: CLEAR,DELETE,SAVE,LOAD, and SAMPLES: DELETE)
function clearStackButton_Callback(hObject, ~, handles) %#ok<DEFNU>

set(handles.cropAboveLabel, 'Enable','on');
set(handles.cropBelowLabel, 'Enable','on');
set(handles.cropAboveLabel,'Enable','on')
set(handles.stackListbox,'String','');
set(handles.pixelWNcalListbox, 'Enable','on');

fields = {'stack','stackIDs','rawDatasetsToPass','selectedDataset'};
for x=1:length(fields)
    try
        handles = rmfield(handles,fields{x});
    catch
    end
end

handles.stack = [];

handles.rawDatasetsToPass = [];
handles.selectedDataset = [];

guidata(hObject, handles);

end
function removeFromStackButton_Callback(hObject, ~, handles) %#ok<DEFNU>

listOfDatasets = get(handles.stackListbox,'String');
selIndex = get(handles.stackListbox,'Value');
selDatasetName = listOfDatasets{selIndex}; % selected dataset name

totalnumSamples = 0; % initialize numSamples
totalnumSamples2 = 0;

for x=1:selIndex
    currSam = str2double(handles.rawDatasetsToPass{x,2});
    totalnumSamples = currSam + totalnumSamples;
end

for x=1:selIndex-1
    currSam = str2double(handles.rawDatasetsToPass{x,2});
    totalnumSamples2 = currSam + totalnumSamples2;
end

% delete the samples from handles.rawDatasetsToPass and handles.stack
handles.rawDatasetsToPass(selIndex,:) = [];

from = totalnumSamples2+1;
to = totalnumSamples;
handles.stack((from:to),:) = [];

% and the stackListbox
listOfDatasets(selIndex) = [];
set(handles.stackListbox,'String',listOfDatasets);
set(handles.stackListbox,'Value',1);

% update feedbackBar
newString = strcat('Deleted dataset! -',{' '},selDatasetName);
set(handles.feedbackBarText,'String',newString);

guidata(hObject, handles); % and save

end
function saveStackButton_Callback(~,~,handles) %#ok<DEFNU>
% check to make sure the stack is not empty before saving
str = get(handles.stackListbox,'String');
if isempty(str)
    return % exit this function if the box is empty
end

% otherwise, save the appropriate handles into a '.mat' file in the
% appropriate directory
dataArrayToPass = cell(1);
S = [];
S.stack = handles.stack;
S.wnn = handles.wnn;
S.rawDatasetsToPass = handles.rawDatasetsToPass;
S.cropAbove = handles.cropAbove;
S.cropBelow = handles.cropBelow;

dataArrayToPass(1) = {S}; %#ok<NASGU>

nameGen = handles.rawDatasetsToPass;
numDatasets = size(nameGen,1);

% Saves files like:
% 'stacked_E-#9-sp9_BSA-2000ug-uL-sp5.mat'
% (sp -  # spectra per dataset)
fileNameGen = '';
for x=1:numDatasets
    if x == 1
        fileNameGen = strcat('stacked_',nameGen(x,1),'-sp',nameGen(x,2));
    elseif x == numDatasets
        fileNameGen = strcat(fileNameGen,'_',nameGen(x,1),'-sp',nameGen(x,2),'.mat');
    else
        fileNameGen = strcat(fileNameGen,'_',nameGen(x,1),'-sp',nameGen(x,2));
    end
end

cd(handles.savedStacksPath);
[file,path] = uiputfile(fileNameGen,'Save file name');
cd(path)
save(file,'dataArrayToPass')
cd(handles.rootPath);

% update feedback bar
newString = strcat('Saved file! -',{' '},fileNameGen);
set(handles.feedbackBarText,'String',newString);
end
function loadStackButton_Callback(hObject, ~, handles) %#ok<DEFNU>
% look for the '.mat' file specified by userinput
cd(handles.savedStacksPath);
[file,path] = uigetfile('*.mat','Load data into stack');
cd(path)
S = load(file);
cd(handles.rootPath);
S = S.dataArrayToPass{1};
handles.stack = S.stack;
handles.wnn = S.wnn;
handles.rawDatasetsToPass = S.rawDatasetsToPass;
handles.cropAbove = S.cropAbove;
handles.cropBelow = S.cropBelow;


guidata(hObject, handles);

newString = strcat('Loaded file! - ',file);
set(handles.feedbackBarText,'String',newString);

set(handles.applyPCAbutton, 'Enable','on');

sampleNames(:,1) = handles.rawDatasetsToPass(:,1);
set(handles.stackListbox,'String',sampleNames);

end
function deleteSample_Callback(hObject, ~, handles) %#ok<DEFNU>
% delete lines from the figures (takes the selections from the
% 'samplesListbox_Callback' function

try
    delete(handles.selectedBGLine);
    delete(handles.selectedRawLine);
catch
    currentSampleNamesInGraph = get(handles.samplesListbox,'String'); % get the list of names
    selectedSampleName = strsplit(currentSampleNamesInGraph{1},{' ','.'}); % get the name of the selected sample
    selectedSampleN = str2double(selectedSampleName{2});
    color = handles.colorMatrix2{selectedSampleN};
    h = findall(0,'type','axes'); % return all of the graphs ('axes') in the figure window
    allLinesBG = findall(h(1),'type','line'); % return the lines in the background graph
    allLinesRawData = findall(h(2),'type','line'); % and raw data graph
    handles.selectedBGLine = findall(allLinesBG, 'Color',color); % and grab the ones matching the color we are looking for
    handles.selectedRawLine = findall(allLinesRawData, 'Color',color);
    delete(handles.selectedBGLine);
    delete(handles.selectedRawLine);
    handles.selectedIndex = 1;
    guidata(hObject, handles);
end

% delete the names of the samples from the listbox
a = get(handles.samplesListbox,'String');
set(handles.samplesListbox,'Value',1);

a(handles.selectedIndex) = [];
set(handles.samplesListbox,'String',a);

% reduce the num of samples text box by 1
numSamples = str2double(get(handles.numSamplesLabel,'String'))-1;
set(handles.numSamplesLabel,'String',numSamples);

%delete the data from the handles.data stream so you can apply new
%backgrounds to the data with deletions
handles.data(handles.selectedIndex,:) = [];
handles.data_mat(handles.selectedIndex,:) = [];
handles.selectedDataset(handles.selectedIndex,:) = [];
handles.sampleNum(handles.selectedIndex,:) = [];
handles = rmfield(handles, 'selectedBGLine');
guidata(hObject, handles);

end

%% analyze Stack
function applyPCAbutton_Callback(hObject, ~, handles) %#ok<DEFNU>

% Run the PCA analysis
[coeff,scores,pcvars] = pca(handles.stack);
handles.coeff = coeff;
handles.scores = scores;
handles.pcvars = pcvars;

% Save handles structure
guidata(hObject, handles);

% and load the first subroutine
StackHQ(handles);

end

%% hanging CreateFcn's
% most of these just set up the white bckgnd for listboxes and popupmenus
function figure1_CreateFcn(~, ~, ~) %#ok<DEFNU>
end
function clearStackButton_CreateFcn(~, ~, ~) %#ok<DEFNU>
end
function numSamplesLabel_CreateFcn(~, ~, ~) %#ok<DEFNU>
end
function rawSpectraAxes_CreateFcn(~, ~, ~) %#ok<DEFNU>
end
% dataset panel:
function datasetListbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function smoothingDSPopupmenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function lagrangeMultEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function pixelWNcalListbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
% background panel:
function plotBGAvgCheckbox_Callback(~, ~, ~) %#ok<DEFNU>

end
function backgroundListbox_CreateFcn(hObject, ~, handles) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
try
    handles.blanks = Func_bg_poly(handles);
catch
end
end
function cropAboveLabel_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function cropBelowLabel_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function smoothingBGPopupmenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function lagrangeMultiEditBG_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function autoFledit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function splitFitAtEditbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
% samples and stacks panel:
function samplesListbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function stackListbox_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function backgroundSpectraEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function rawSpectraAxesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%% TO IMPLEMENT
function permanentlyDeleteSampleButton_Callback(~, ~, ~) %#ok<DEFNU>
end
function upDirButton_Callback(~, ~, ~) %#ok<DEFNU>
end


function backgroundSpectraSaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.backgroundSpectraEdit,'String'),'.pdf'];
 export_fig(handles.backgroundSpectraAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function backgroundSpectraEdit_Callback(hObject, eventdata, handles)
end
function rawSpectraAxesSaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.rawSpectraAxesEdit,'String'),'.pdf'];
 export_fig(handles.rawSpectraAxes, fileName,'-q101','-rgb','-painters','-dpdf');

end
function rawSpectraAxesEdit_Callback(hObject, eventdata, handles)
end
function normalizeDataCheckbox_Callback(hObject, eventdata, handles)
end
