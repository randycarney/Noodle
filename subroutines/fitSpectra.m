% Begin initialization code - DO NOT EDIT
function varargout = fitSpectra(varargin)
% Edit the above text to modify the response to help fitSpectra

% Last Modified by GUIDE v2.5 20-Oct-2016 14:34:46
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @fitSpectra_OpeningFcn, ...
    'gui_OutputFcn',  @fitSpectra_OutputFcn, ...
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

%% pre-load
function varargout = fitSpectra_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end
function fitSpectra_OpeningFcn(hObject, ~, handles, varargin)

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

passedInHandles = varargin{1};
handles.stack = passedInHandles.stack;
handles.wnn = passedInHandles.wnn;
handles.rawDatasets = passedInHandles.rawDatasets;
handles.colorMatrix = passedInHandles.colorMatrix;
handles.colorMatrix2 = passedInHandles.colorMatrix2;
handles.colorMatrixRGBtrip = passedInHandles.colorMatrixRGBtrip;
handles.coeff = passedInHandles.coeff;
handles.scores = passedInHandles.scores;
handles.pcvars = passedInHandles.pcvars;
handles.cropAbove = passedInHandles.cropAbove;
handles.cropBelow = passedInHandles.cropBelow;
handles.clusterNames = passedInHandles.clusterNames;
handles.numClusters = passedInHandles.numClusters;
handles.stax = passedInHandles.stax;
handles.rootPath = passedInHandles.rootPath;
handles.standardsPath = passedInHandles.standardsPath;

guidata(hObject,handles); % and then save them

set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

% Load all the standard files into the listbox
cd(handles.standardsPath); % change directory to the rawData path
standardNames = dir('*.SPE'); % get the file info for every .mat file
standardNames = {standardNames(~[standardNames.isdir]).name}; % and extract the string names
rawStandardNames = standardNames'; % rotate the matrix bc I think better in columns for some reason
if length(rawStandardNames) < 1
    rawStandardNames = '(no standards on file)';
end
handles.rawStandardNames = rawStandardNames;
guidata(hObject,handles); % and then save them
set(handles.standardsMenu,'String',rawStandardNames)
set(handles.clusterMenu,'String',handles.clusterNames) % set up clusterMenu1
set(handles.clusterMenu,'Value',1)
set(handles.standardsMenu,'Value',5)
handles.numClusters = 1;
guidata(hObject,handles); % and then save them

standardsMenu_Callback(handles.standardsMenu, 1, handles)

axes(handles.clusterAxes);
cla;
set(gca,'Tag','clusterAxes')
ylabel('Intensity (A.U.)');
xlabel('Rel. wavenumber (cm^-^1)')
hold on
hline(0)
for x=1:1
    tf = get(handles.overlayCheckbox,'Value');
    if tf == 1
        shadedErrorBar(handles.wnn,mean(handles.stack,1),std(handles.stack,1),'k',2)
        plot(handles.wnn,mean(handles.stack,1),'k','LineWidth',1,'LineStyle','--')
    end
    plot(handles.wnn,(mean(handles.stax(x).stack,1)-mean(handles.stack,1)),handles.colorMatrix{x},'LineWidth',2)
    set(gca,'xlim',[handles.cropBelow handles.cropAbove])
end

cd(handles.rootPath)
end

function clusterMenu_Callback(hObject, ~, handles) %#ok<DEFNU>

contents = cellstr(get(hObject,'String')); %returns testy contents as cell array
selCluster = strsplit(contents{get(hObject,'Value')},{'Cluster ','</b>'}); %returns selected item from testy
numClusters = str2double(selCluster{2});
handles.numClusters = numClusters;
guidata(hObject,handles); % and then save them
plotNewCluster(1,numClusters,handles)
standardsMenu_Callback(handles.standardsMenu, 1, handles)

end

function plotNewCluster(~,n,handles)

axes(handles.clusterAxes)
cla
set(gca,'Tag','clusterAxes')
hold on
tf = get(handles.overlayCheckbox,'Value');
if tf == 1
    shadedErrorBar(handles.wnn,mean(handles.stack,1),std(handles.stack,1),'k',2)
    plot(handles.wnn,mean(handles.stack,1),'k','LineWidth',1,'LineStyle','--')
end
hline(0)
plot(handles.wnn,(mean(handles.stax(n).stack,1)-mean(handles.stack,1)),handles.colorMatrix{n},'LineWidth',2)
set(gca,'xlim',[handles.cropBelow handles.cropAbove])

end

function overlayCheckbox_Callback(~, ~, handles) %#ok<DEFNU>
if any(strcmp('numClusters',fieldnames(handles))) == 1
    numClusters = handles.numClusters;
else
    numClusters = 1;
end
plotNewCluster(1,numClusters,handles)
end

function standardsMenu_Callback(hObject, ~, handles)

cd(handles.standardsPath);

selectedStandards = get(handles.standardsMenu,'Value');

if selectedStandards == 0
    selectedStandards = 1;
end

sel = length(selectedStandards);
handles.selectedStandards = selectedStandards;
guidata(hObject, handles);

% make a selected standards name for the figure legend
allStandards = get(handles.standardsMenu,'String');
for row=1:sel
    currName = allStandards{handles.selectedStandards(row)};
    % corrName = strsplit(currName,' ');
    selectedNames{row} = currName;
end

% 1. get the user input crop values for wavenumber (wn -> wnn)
cropAbove = handles.cropAbove;
cropBelow = handles.cropBelow;
guidata(hObject, handles); % and save

wnn = handles.wnn;

for N = 1:sel
    % 1. fill a matrix with the background names
    % (currently unused at the moment)
    selectedBgNames(1,N) = cellstr(handles.rawStandardNames{selectedStandards(1,N)}); % create a cell array of bg names
    % 2. read in the current background SPE file in blankPlaceholder
    standardsSpectra = double(readSPE(handles.rawStandardNames{selectedStandards(1,N)})); % read each bg SPE into a matrix
    % 3. manually 'accumulate' the scans if the user did not specificy
    % accumulations during aquisition
    if size(standardsSpectra,3) > 1
        standardsSpectra = sum(squeeze(standardsSpectra),2)';
    end
    % blankSpectra should now be of dimension: [1 x 1340]
    % 4. It is cropped according to the user-input values (passed in
    % handles)
    croppedBlankSpectra = standardsSpectra(wnn>=cropBelow & wnn<=cropAbove);
    % 5. And finally smoothed
    % grab the Lagrange parameter from the user input box
    bgLG = str2double(get(handles.lagrangeStandardsEdit,'String'));
    qz1 = whittaker_smoother(croppedBlankSpectra,bgLG);
    
    % and appended to a loopholder matrix that gets passed back as this
    % function's output
    blankMat(N,:) = qz1;    % append the smoothed, chopped blank to a matrix
    
end

try
    handles.blankMat = blankMat;
    guidata(hObject, handles);
catch
    cd(handles.rootPath)
    return
end
% finish by appending a n-th order% polynomial function representing
% system autofluorescence (will be used to fit to the data with asls)
% polynomial order is specified by user input:
autoFl = str2double(get(handles.autoFlEdit,'String'));

% this matrix is output to the return of the called function
blankMat2 = [blankMat;mypoly(autoFl,length(blankMat))];

axes(handles.standardsAxes)
cla
set(gca,'Tag','standardsAxes');
hold on
set(gca, 'xlim',[cropBelow cropAbove]);
ylim auto
if sel > 1
    for x=1:sel
        standardLines(x) = plot(wnn,blankMat(x,:),'Color',handles.colorMatrix2{x});
    end
    legend(standardLines,selectedNames,'FontSize',10)
else
    standardLines = plot(wnn,blankMat,'Color',rgb('Green'));
    legend(standardLines,selectedNames,'FontSize',10)
end

axes(handles.standardsFitAxes)
cla; hold on;
set(gca,'Tag','standardsFitAxes')
set(gca, 'xlim',[cropBelow cropAbove]);

val = get(handles.clusterMenu,'Value');
test1 = (mean(handles.stax(val).stack,1)-mean(handles.stack,1));

% include or don't include the autofl. poly in blanks (disableAFcheckbox)
if get(handles.disableAFcheckbox,'Value') ~= 1 
    [x,~,~,~] = lscov(blankMat2',test1');
    fit = (x'*blankMat2);
    corrected = test1-(x'*blankMat2);
    if get(handles.plotPoly,'Value') ~= 1
        % both boxes unchecked
        for N=1:sel
            toPlot = x(N)*blankMat(N,:);
            lines(N) = plot(wnn,toPlot,'Color',handles.colorMatrix2{N},'LineStyle','--');
            formatSpec = [selectedNames{N},', fit: %+5.4f'];
            newStr = sprintf(formatSpec,x(N));
            selectedNames(N) = cellstr(newStr);
        end        
        legend(lines,selectedNames,'FontSize',11)
    else
        % poly autofl. fit enabled, plotting them disabled
        for N=1:size(blankMat2,1)
            toPlot = x(N)*blankMat2(N,:);
            lines(N) = plot(wnn,toPlot,'Color',handles.colorMatrix2{N},'LineStyle','--');
        end
        lines(sel+1:end) = [];
        legend(lines,selectedNames,'FontSize',11)
    end
else
    % both boxes checked
    [x,~,~,~] = lscov(blankMat',test1');
    fit = (x'*blankMat);
    corrected = test1-(x'*blankMat);
    standardLines = plot(wnn,fit,'Color',rgb('Green'),'LineStyle','--');
    legend(standardLines,selectedNames,'FontSize',11)
end
hline(0)
%
% coeffs = asls(blankMat2',test1',.05,[]); % estimate the background coefficients using AsLS
% corrected = test1-(coeffs'*blankMat2);  % construct the estimated background and subtract it from data

axes(handles.fitAxes)
cla; hold on;
set(gca,'Tag','fitAxes')
set(gca, 'xlim',[cropBelow cropAbove]);
L1 = plot(wnn,test1,'Color',handles.colorMatrix{handles.numClusters},'LineWidth',1);
L2 = plot(wnn,fit,'Color','k','LineWidth',1,'LineStyle','-.');
hline(0)
legend([L1;L2],{'cluster avg. diff.';'best fit to cluster'},'FontSize',11)

axes(handles.residualsAxes)
cla; hold on;
set(gca,'Tag','residualsAxes')
set(gca, 'xlim',[cropBelow cropAbove]);
L1 = plot(wnn,corrected,'Color','k');
xlabel('Rel. wavenumber (cm^-^1)')
hline(0)
legend(L1,{'best fit residuals'},'FontSize',11)

cd(handles.rootPath);

end

function autoFlEdit_Callback(hObject, eventdata, handles) %#ok<DEFNU>
standardsMenu_Callback(hObject, eventdata, handles)
end
function lagrangeStandardsEdit_Callback(~, ~, ~) %#ok<DEFNU>
end
function disableAFcheckbox_Callback(hObject, eventdata, handles) %#ok<DEFNU>

standardsMenu_Callback(hObject, eventdata, handles)
end
function plotPoly_Callback(hObject, eventdata, handles) %#ok<DEFNU>
if get(hObject,'Value') == 1
    set(handles.disableAFcheckbox,'Value',0);
end
standardsMenu_Callback(hObject, eventdata, handles)
end

%% hanging createFcns
function lagrangeStandardsEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function autoFlEdit_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function standardsMenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function clusterMenu_CreateFcn(hObject, ~, ~) %#ok<DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function scaledStandardsEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function standardsFitOverlayEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function fitResidualsEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function standardsSpectraEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ScaledStandardsSaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.scaledStandardsEdit,'String'),'.pdf'];
 export_fig(handles.standardsFitAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function scaledStandardsEdit_Callback(hObject, eventdata, handles)
end

function standardsFitOverlaySaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.standardsFitOverlayEdit,'String'),'.pdf'];
 export_fig(handles.fitAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function standardsFitOverlayEdit_Callback(hObject, eventdata, handles)
end

function fitResidualsSaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.fitResidualsEdit,'String'),'.pdf'];
 export_fig(handles.residualsAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function fitResidualsEdit_Callback(hObject, eventdata, handles)
end

function standardsSpectraSaveButton_Callback(hObject, eventdata, handles)
fileName = [get(handles.standardsSpectraEdit,'String'),'.pdf'];
 export_fig(handles.standardsAxes, fileName,'-q101','-rgb','-painters','-dpdf');
end
function standardsSpectraEdit_Callback(hObject, eventdata, handles)
end
