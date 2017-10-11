% Begin initialization code - DO NOT EDIT
function varargout = PCanalysis(varargin)
clc
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PCanalysis_OpeningFcn, ...
    'gui_OutputFcn',  @PCanalysis_OutputFcn, ...
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

% --- Outputs from this function are returned to the command line.
function varargout = PCanalysis_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

% --- Executes just before PCanalysis is made visible.
function PCanalysis_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for samplePicker3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

%% VARARGIN is used to pass in the raw spectra data generated from the 'samplePicker' window
% 1. stack: the raw dataset spectra
% 2. wnn: the user-clipped wavenumber distribution
% 3. rawDatasets: the names, numbers, and colors of datasets

passedInHandles = varargin{1};
handles.stack = passedInHandles.stack;
handles.wnn = passedInHandles.wnn;
handles.rawDatasets = passedInHandles.rawDatasetsPassed;
handles.colorMatrix = passedInHandles.colorMatrix;
handles.colorMatrix2 = passedInHandles.colorMatrix2;
handles.colorMatrixRGBtrip = passedInHandles.colorMatrixRGBtrip;
handles.coeff = passedInHandles.coeff;
scores = passedInHandles.scores;

handles.scores = passedInHandles.scores;
handles.pcvars = passedInHandles.pcvars;

handles.cropAbove = passedInHandles.cropAbove;
handles.cropBelow = passedInHandles.cropBelow;

guidata(hObject,handles); % and then save them


%% set up the colors
handles.colorMatrix = {'g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y','g','r','b','k','m','c','y'};
colorMatrix2 = {'Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen','Green','Red','Blue','Black','Purple','Turquoise','Gold','DarkRed','Gray','LightSlateGray','Brown','OrangeRed','DeepSkyBlue','Indigo','SeaGreen'};
for x=1:length(colorMatrix2)
    colorMatrix2{x} = rgb(colorMatrix2(x));
end
handles.colorMatrix2 = colorMatrix2;
handles.colorMatrixRGBtrip = {[0 1 0],[1 0 0],[0 0 1],[0 0 0],[1 0 1],[0 1 1],[1 1 0],[1 0 0],[0 1 0],[0 0 1],[0 0 0],[1 0 1],[0 1 1],[1 1 0]};

guidata(hObject,handles); % and then save them

%% set up the drop box boxes that list each PC
numOfScores = size(scores,1)-1;
PCcatGen = @(number) ['PC',number,''];

for x=1:numOfScores
    listOfScores{x} = PCcatGen(num2str(x));
end

set(handles.PC1popupmenu,'String',listOfScores);
set(handles.PC1popupmenu,'Value',1);
set(handles.PC2popupmenu,'String',listOfScores);
set(handles.PC2popupmenu,'Value',2);
set(handles.PC3popupmenu,'String',listOfScores);
set(handles.PC3popupmenu,'Value',3);

%% Plot the average spectra and first 4 PC coeffs
for n=1:4
    if n == 4
        axes(handles.globalAverage)
        cla
        set(gca,'Tag','globalAverage')
        shadedErrorBar(handles.wnn,mean(handles.stack,1),std(handles.stack,1),'k',2)
        hold on
        plot(handles.wnn,mean(handles.stack,1),'k','LineWidth',1)
        plotStandards()
        set(gca,'xlim',[handles.cropBelow handles.cropAbove])
        ylabel('Intensity (A.U.)');
        xlabel('Rel. wavenumber (cm^-^1)')
    else
        field = ['PC',num2str(n),'coeffs'];
        axes(handles.(field))
        cla
        set(gca,'Tag',field)
        plotCoeffs(n)
        plotStandards()
        set(gca,'xlim',[handles.cropBelow handles.cropAbove])
    end
end

    function plotCoeffs(x)
        plot(handles.wnn,handles.coeff(:,x),handles.colorMatrix{x})
        ylabel('Intensity (A.U.)');
        xlabel('Rel. wavenumber (cm^-^1)')
        hold on
        hline(0)
    end

%% Plot PC scores

colergen1 = @(color1,color2,color3,text) ['<html><table border=0 width=400 font color=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];
colergen2 = @(color1,color2,color3,text) ['<html><table border=0 width=400 bgcolor=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];

colorsss = cellfun(@(x) x*255,colorMatrix2','un',0); % convert to rgb style for html

for l=1:size(handles.rawDatasets,1)
    data{l,1} = colergen2(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),'');
    data{l,2} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasets{l,2}); % set the #samples column
    data{l,3} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasets{l,1}); % set the dataset column
end

set(handles.samplesTable, 'Data',data)

dummyYvaluesOnes = ones(size(scores,1), 1);
dummyYvalues = dummyYvaluesOnes;

textSpacing = 0.1;
dummyYvaluesText = [];
% make the numbers and color distributions for the input datasets
sampleNumbers = [];
for x=1:size(handles.rawDatasets,1)
    numOfSamplesInThisDataset = str2double(handles.rawDatasets{x,2});
    numOfSamples(x) = numOfSamplesInThisDataset;
    currentLength = size(sampleNumbers,1);
    for y=1:numOfSamplesInThisDataset
        sampleNumbers(currentLength+y,:) = y;
        colorMatrixFromData(currentLength+y,:) = handles.colorMatrix2{x};
        dummyYvaluesText(currentLength+y,:) = (1.05+(x*textSpacing))*dummyYvaluesOnes(currentLength+y,:);
    end
end
size(colorMatrixFromData,1)

sampleNumbers = num2str(sampleNumbers);
numOfSamples = sum(numOfSamples); % then sum them up

handles.sampleNumbers = sampleNumbers;
handles.numOfSamples = numOfSamples;
handles.colorMatrixFromData = colorMatrixFromData;
handles.dummyYvaluesText = dummyYvaluesText;
handles.dummyYvalues = dummyYvalues;
guidata(hObject, handles);

axes(handles.PC1scores)
cla
plotScatter(1)

axes(handles.PC2scores)
cla
plotScatter(2)

axes(handles.PC3scores)
cla
plotScatter(3)

    function plotScatter(x)
        boxplot(scores(:,x),'orientation','horizontal')
        hold on
        scatter(scores(:,x),dummyYvalues,30,colorMatrixFromData, 'filled')
        for z=1:length(colorMatrixFromData)
            text(scores(z,x),dummyYvaluesText(z,:),sampleNumbers(z,:),'Color',colorMatrixFromData(z,:),'FontSize',10)
        end
        hold on
        vline(0)
        set(gca, 'YTickLabel',[]);
        set(gca, 'YLim',[0.9 2]);
        xlabel('Rel. wavenumber (cm^-^1)')
    end

handles.didScaleTogetherBool = 0;
guidata(hObject, handles);

end

function PC1popupmenu_Callback(hObject, ~, handles)
redrawPlot(hObject,handles,1)
end
function PC2popupmenu_Callback(hObject, ~, handles)
redrawPlot(hObject,handles,2)
end
function PC3popupmenu_Callback(hObject, ~, handles)
redrawPlot(hObject,handles,3)
end

function redrawPlot(hObject,handles,selectedGraph)

switch selectedGraph
    case 1
        axes(handles.PC1coeffs)
    case 2
        axes(handles.PC2coeffs)
    case 3
        axes(handles.PC3coeffs)
end

currentListOfScores = cellstr(get(hObject,'String'));
currentSelection = get(hObject,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
x = str2double(split{2});
cla
% plot coeff distribution
plot(handles.wnn,handles.coeff(:,x),handles.colorMatrix{x})
hold on

plotStandards()

switch selectedGraph
    case 1
        axes(handles.PC1scores)
    case 2
        axes(handles.PC2scores)
    case 3
        axes(handles.PC3scores)
end

cla
scatter(handles.scores(:,x),handles.dummyYvalues,30,handles.colorMatrixFromData, 'filled')
for z=1:length(handles.colorMatrixFromData)
    text(handles.scores(z,x),handles.dummyYvaluesText(z,:),handles.sampleNumbers(z,:),'Color',handles.colorMatrixFromData(z,:),'FontSize',10)
end
hold on
vline(0)
boxplot(handles.scores(:,x),'orientation','horizontal')
set(gca, 'YTickLabel',[]);
set(gca, 'YLim',[0.9 2]);
hold on

end

%% scale buttons
function scaleTogetherButton_Callback(hObject, ~, handles)

if handles.didScaleTogetherBool == 0
    handles.XLim1 = get(handles.PC1scores, 'XLim');
    handles.XLim2 = get(handles.PC2scores, 'XLim');
    handles.XLim3 = get(handles.PC3scores, 'XLim');
    
    
    globalXmin = min([handles.XLim1(:); handles.XLim2(:); handles.XLim3(:)]);
    globalXmax = max([handles.XLim1(:); handles.XLim2(:); handles.XLim3(:)]);
    XLim = [globalXmin globalXmax];
    
    set(handles.PC1scores, 'XLim',XLim);
    set(handles.PC2scores, 'XLim',XLim);
    set(handles.PC3scores, 'XLim',XLim);
    
    handles.didScaleTogetherBool = 1;
    guidata(hObject,handles);
end

end
function scaleIndepButton_Callback(hObject, ~, handles) %#ok<*DEFNU>

if handles.didScaleTogetherBool == 1
    set(handles.PC1scores, 'XLim',handles.XLim1);
    set(handles.PC2scores, 'XLim',handles.XLim2);
    set(handles.PC3scores, 'XLim',handles.XLim3);
    handles.didScaleTogetherBool = 0;
    guidata(hObject,handles);
end

end

%% hanging CreateFcn's
% most of these just set up the white bckgnd for listboxes and popupmenus
function listbox1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function PC1popupmenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function PC2popupmenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function PC3popupmenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function PCscoresEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function PCcoeffEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function globalAvgEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%% To Implement

function plotStandards()
% this should highlight areas for standards - currently not implemented
% Ylim = get(gca,'ylim');
% bottomY = Ylim(1);
% topY = Ylim(2);
% heightY = Ylim(2)-bottomY;
% % plot protein
% 
% protein = [820,30; 950,70; 1200,170; 1435,65; 1640,50];
% lipid = [1250,100; 1440,35; 1640,25];
% nucleicAcid = [800,140; 1025,125; 1275,165; 1540,60];

% for z=1:size(protein,1)
%     rectangle('Position',[protein(z,1),bottomY,protein(z,2),(topY-heightY)],'FaceColor',rgb('LightCyan'),'EdgeColor','none')
% end
% for z=1:size(lipid,1)
%     rectangle('Position',[lipid(z,1),(1/3)*topY,lipid(z,2),heightY],'FaceColor',rgb('GreenYellow'),'EdgeColor','none')
% end
% for z=1:size(nucleicAcid,1)
%     rectangle('Position',[nucleicAcid(z,1),(2/3)*topY,nucleicAcid(z,2),heightY],'FaceColor',rgb('Honeydew'),'EdgeColor','none')
% end

end

function globalAvgSaveButton_Callback(~, ~, handles)
fileName = [get(handles.globalAvgEdit,'String'),'.pdf'];
 export_fig(handles.globalAverage, fileName,'-q101','-rgb','-painters','-dpdf');
end
function PCcoeffSaveButton_Callback(~, ~, handles)
fileName = [get(handles.PCcoeffEdit,'String'),'.pdf'];
 export_fig(handles.PC1coeffs, fileName,'-q101','-rgb','-painters','-dpdf');
end
function PCscoresSaveButton_Callback(~, ~, handles)
fileName = [get(handles.PCscoresEdit,'String'),'.pdf'];
 export_fig(handles.PC1scores, fileName,'-q101','-rgb','-painters','-dpdf');
end
function PCcoeffEdit_Callback(~, ~, ~)
end
function PCscoresEdit_Callback(~, ~, ~)
end
function globalAvgEdit_Callback(~, ~, ~)
end

