% Begin initialization code - DO NOT EDIT
function varargout = StackHQ(varargin)
% STACKHQ MATLAB code for StackHQ.fig
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @StackHQ_OpeningFcn, ...
    'gui_OutputFcn',  @StackHQ_OutputFcn, ...
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

function varargout = StackHQ_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

%% --- Executes just before StackHQ is made visible.
function StackHQ_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args

% Choose default command line output for samplePicker3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

%% VARARGIN is used to pass in the raw spectra data generated from the 'samplePicker' window

passedInHandles = varargin{1};
handles.stack = passedInHandles.stack;
handles.wnn = passedInHandles.wnn;
handles.rawDatasetsPassed = passedInHandles.rawDatasetsToPass;
handles.colorMatrix = passedInHandles.colorMatrix;
handles.colorMatrix2 = passedInHandles.colorMatrix2;
handles.colorMatrixRGBtrip = passedInHandles.colorMatrixRGBtrip;
handles.coeff = passedInHandles.coeff;
handles.scores = passedInHandles.scores;
handles.rootPath = passedInHandles.rootPath;
handles.standardsPath = passedInHandles.standardsPath;
handles.pcvars = passedInHandles.pcvars;
handles.cropAbove = passedInHandles.cropAbove;
handles.cropBelow = passedInHandles.cropBelow;

guidata(hObject,handles); % and then save them

% housekeeping
% close figure 1 if it's open
if(any(findall(0,'Type','Figure')==1))
    close(1);
end

%% populate up the datasetTable
updateDatasetTable(handles)
updateSamplesTable(1,handles)

%% plot the total sample averages

% 1. Contruct the averages plot
range = plotAverages(handles);
handles.range = range; % save in the range table
guidata(hObject,handles);

% 2. Plot the first dataset by default in the upper graph
plotSpectra(1,handles)

% set the # of coeffs
set(handles.PCAscoresText,'String',size(handles.coeff,2))

end

function updateDatasetTable(handles)
%1. set colors for each dataset
% (need to use HTML here, it's the only way to color table cells like I want)
colergen1 = @(color1,color2,color3,text) ['<html><table border=0 width=400 font color=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];
colergen2 = @(color1,color2,color3,text) ['<html><table border=0 width=400 bgcolor=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];

colorsss = cellfun(@(x) x*255,handles.colorMatrix2','un',0); % convert to rgb style for html

numDatasets = size(handles.rawDatasetsPassed,1);

data = cell(numDatasets,3);

for l=1:numDatasets
    data{l,1} = colergen2(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),'');
    data{l,2} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasetsPassed{l,2}); % set the #samples column
    data{l,3} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasetsPassed{l,1}); % set the dataset column
end

set(handles.datasetTable,'Data',data)
set(handles.deleteSelDataset,'Enable','off')

end

function updateSamplesTable(n,handles)

try
    datasetName = handles.rawDatasetsPassed{n};
catch
    data{1,2} = 'no samples in dataset';
    set(handles.sampleTable,'Data',data)
    return
end

%1. set colors for each dataset
% (need to use HTML here, it's the only way to color table cells like I want)
colergen1 = @(color1,color2,color3,text) ['<html><table border=0 width=400 font color=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];
colergen2 = @(color1,color2,color3,text) ['<html><table border=0 width=400 bgcolor=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];
colorsss = cellfun(@(x) x*255,handles.colorMatrix2','un',0); % convert to rgb style for html

numSamples = str2double(handles.rawDatasetsPassed{n,2});
data = cell(numSamples,2);

for l=1:numSamples
    nameGen = {datasetName,num2str(l)};
    name = strjoin(nameGen);
    data{l,2} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),name); % set the dataset column
    data{l,1} = colergen2(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),'');
end

set(handles.sampleTable,'Data',data)
set(handles.deleteSampleButton,'Enable','off')

end

function rangeToUse = plotAverages(handles)

% regenerates the globalMeanAxes using:
% 1. handles.rawDatasetsPassed
% 2. handles.range
% 3. handles.stack

numDatasets = size(handles.rawDatasetsPassed,1);

totalNumSamples = 0;
rangeToUse = zeros(numDatasets,2);
for x=1:numDatasets
    rangeToUse(x,1) = totalNumSamples+1;
    numSamples = str2double(handles.rawDatasetsPassed{x,2});
    totalNumSamples = totalNumSamples + numSamples;
    rangeToUse(x,2) = totalNumSamples;
end

axes(handles.globalMeanAxes)
cla
set(gca,'Tag','globalMeanAxes') % cla also clears the Tag, so need to re-state it

% iterate through each dataset and calculate/plot the mean + s.d.
for x=1:size(rangeToUse,1)
    from = rangeToUse(x,1);
    to = rangeToUse(x,2);
    numSamples = to-(from-1);
    selSampleSpectra(1:numSamples,:) = handles.stack(from:to,:);
    data_mean = mean(selSampleSpectra,1);
    data_std = std(selSampleSpectra);
    if length(handles.wnn) == length(data_std)
        shadedErrorBar(handles.wnn,data_mean,data_std,handles.colorMatrix2{x},2)
    else
    end
    hold on
end
ylabel('Intensity (A.U.)');
xlabel('Rel. wavenumber (cm^-^1)')
set(gca,'xlim',[handles.cropBelow handles.cropAbove])

end

function plotSpectra(n,handles)

try
    datasetName = handles.rawDatasetsPassed{n};
catch
    axes(handles.selDatasetSpectraAxes)
    cla
    set(handles.datasetName,'String','no datasets in Stack')
    return
end

set(handles.datasetName,'String',datasetName)
set(handles.datasetName,'ForegroundColor',handles.colorMatrix2{n})

axes(handles.selDatasetSpectraAxes)
cla
hold on
set(gca,'Tag','selDatasetSpectraAxes') % cla also clears the Tag, so need to re-state it

range = handles.range(n,:);
from = range(1);
to = range(2);

for x=from:to
    plot(handles.wnn, handles.stack(x,:),'Color', handles.colorMatrix2{(x-(from-1))})
end
ylabel('Intensity (A.U.)');
xlabel('Rel. wavenumber (cm^-^1)')
set(gca,'xlim',[handles.cropBelow handles.cropAbove])

end

function datasetTable_CellSelectionCallback(hObject, eventdata, handles) %#ok<*DEFNU>

index = eventdata.Indices;

if isempty(index)
    n = 1;
    handles.selDatasetRow = n; % save globally
    guidata(hObject,handles); % and then save them
    plotSpectra(n,handles)
    set(handles.deleteSelDataset,'Enable','off')
else
    n = index(1);
    handles.selDatasetRow = n; % save globally
    guidata(hObject,handles); % and then save them
    plotSpectra(n,handles)
    set(handles.deleteSelDataset,'Enable','on')
    updateSamplesTable(n,handles)
end

end

function sampleTable_CellSelectionCallback(hObject, eventdata, handles)

index = eventdata.Indices;

if ~isempty(index)
    n = index(1);
else
    n = 1;
    set(handles.deleteSampleButton,'Enable','off')
end
handles.selSampleRow = n; % save globally
guidata(hObject,handles); % and then save them

% b. get the color that this sample would have been graphed in
color = handles.colorMatrix2{n};

% c. access the figure and grab the appropriate colored line in graph
h = findobj(0,'Tag','selDatasetSpectraAxes');
% first remove any bolded lines from previous selections
allLines = findall(h,'type','line'); % return all graphed lines in both graphs
set(allLines,'LineWidth',0.5); % set all of ther widths to not-bold (0.5)
% then find the lines we want to bold
allLinesRawData = findall(h,'type','line'); % and raw data graph
selectedRawLine = findall(allLinesRawData, 'Color',color);
% Bold the selected lines
set(selectedRawLine,'LineWidth',2);

% 3. Finally store the selected lines in the GUIDATA handles for deletion if
% necessary
handles.selectedRawLine = selectedRawLine;
guidata(hObject, handles);

set(handles.deleteSampleButton,'Enable','on')

end

%% DATA HANDLING (DATASET: DELETE)
function deleteSelDataset_Callback(hObject, ~, handles)

% close figure 1 if it's open
if(any(findall(0,'Type','Figure')==1)),
    close(1);
    wasOpen = 1;
else
    wasOpen = 0;
end;

% upon dataset deletion, need to update the:
% 1. handles.rawDatasetsPassed
% 2. handles.range
% 3. handles.stack

% get the selected Row
try
    r = handles.selDatasetRow;
catch
    return
end

% and then delete the data appropriate handles
from = handles.range(r,1);
to = handles.range(r,2);

handles.range(r,:) = [];
handles.rawDatasetsPassed(r,:) = [];
handles.stack(from:to,:) = [];

guidata(hObject,handles); % and then save them

% finally re-draw everything and re-do the PCA

% 1. re-draw everything
%% populate up the datasetTable
updateDatasetTable(handles)
updateSamplesTable(r,handles)

%% plot the total sample averages

% Contruct the averages plot
range = plotAverages(handles);
handles.range = range; % save in the range table
guidata(hObject,handles);

% Then plot the first dataset by default in the upper graph
plotSpectra(1,handles)

% 2. re-do PCA

tf = get(handles.varNormCheckbox,'Value'); % get state of 'variance normalized' checkbox
if tf == 1
    [coeff,scores,pcvars] = pca(zscore(handles.stack));
else
    [coeff,scores,pcvars] = pca(handles.stack);
end
handles.coeff = coeff;
handles.scores = scores;
handles.pcvars = pcvars;
set(handles.PCAscoresText,'String',size(coeff,2))
guidata(hObject, handles); % and save it

if wasOpen ==1
    figure(1)
    title('All 1D PCA boxplots')
    hold on
    boxplot(handles.scores,'orientation','horizontal')
    hold on
end

end

function deleteSampleButton_Callback(hObject, eventdata, handles)

% close figure 1 if it's open
if(any(findall(0,'Type','Figure')==1)),
    close(1);
    wasOpen = 1;
else
    wasOpen = 0;
end;

% 1. delete the data from the stack
numSamplesSoFar = 0;
for x=1:handles.selDatasetRow
    if x == handles.selDatasetRow
         finalRow = handles.selSampleRow + numSamplesSoFar;
         size(handles.stack,1)
         handles.stack(finalRow,:) = [];
         size(handles.stack,1)
         guidata(hObject,handles);
    else
        numSamplesSoFar = numSamplesSoFar + str2double(handles.rawDatasetsPassed{x,2});
    end
end

% Run the PCA analysis
tf = get(handles.varNormCheckbox,'Value'); % get state of 'variance normalized' checkbox
if tf == 1
    [coeff,scores,pcvars] = pca(zscore(handles.stack));
else
    [coeff,scores,pcvars] = pca(handles.stack);
end
handles.coeff = coeff;
handles.scores = scores;
handles.pcvars = pcvars;
guidata(hObject,handles);
% set the # of coeffs
set(handles.PCAscoresText,'String',size(handles.coeff,2))

% 2. update the rawDatasetsToPass numSamples column
clc
handles.rawDatasetsPassed(handles.selDatasetRow,2) = cellstr(num2str(str2double(handles.rawDatasetsPassed{handles.selDatasetRow,2})-1));
guidata(hObject,handles);

if (str2double(handles.rawDatasetsPassed{handles.selDatasetRow,2})-1) == 0
    deleteSelDataset_Callback(hObject, eventdata, handles)
end

% 3. Replot the tables/figures
updateDatasetTable(handles)
n = handles.selDatasetRow;
updateSamplesTable(n,handles)

range = plotAverages(handles);
handles.range = range; % save in the range table
guidata(hObject,handles);

plotSpectra(n,handles)

% 4. delete the actual line in the graph

delete(handles.selectedRawLine);

% re-plot dendrogram if necessary
if wasOpen ==1
    figure(1)
    title('All 1D PCA boxplots')
    hold on
    boxplot(handles.scores,'orientation','horizontal')
    hold on
end

end

function varNormCheckbox_Callback(hObject, ~, handles)

tf = get(hObject,'Value'); % get state of 'variance normalized' checkbox
if tf == 1
    [coeff,scores,pcvars] = pca(zscore(handles.stack));
else
    [coeff,scores,pcvars,~,explained,~] = pca(handles.stack);
    disp(explained)
end
handles.coeff = coeff;
handles.scores = scores;
handles.pcvars = pcvars;

guidata(hObject,handles);

set(handles.PCAscoresText,'String',size(coeff,2))

end
function ExportAvgDataFigButton_Callback(~, ~, handles)
fileName = [get(handles.datasetAvgSpecEdit,'String'),'.pdf'];
 export_fig(handles.globalMeanAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function datasetAvgSpecEdit_Callback(~, ~, ~)
end
function AllSpectraExportFigButton_Callback(~, ~, handles)
fileName = [get(handles.AllSpecEdit,'String'),'.pdf'];
 export_fig(handles.selDatasetSpectraAxes,fileName,'-q101','-rgb','-painters','-dpdf');
end
function AllSpecEdit_Callback(~, ~, ~)
end

%% OPEN NEW WINDOWS
function allPCsBoxplotButton_Callback(~, ~, handles)
% plot the box plot in figure(1)
figure(1)
title('All 1D PCA boxplots')
hold on
boxplot(handles.scores,'orientation','horizontal')
hold off
end
function applyPCA2Dbutton_Callback(~, ~, handles)
PCvsPCgen(handles);
end
function applyPCA1Dbutton_Callback(~, ~, handles)
PCanalysis(handles);
end

%% hanging CreateFcn's
% most of these just set up the white bckgnd for listboxes and popupmenus
function PCAscoresText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function datasetAvgSpecEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function AllSpecEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
