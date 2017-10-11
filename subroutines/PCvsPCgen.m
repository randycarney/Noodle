% Begin initialization code - DO NOT EDIT
function varargout = PCvsPCgen(varargin)

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @PCvsPCgen_OpeningFcn, ...
    'gui_OutputFcn',  @PCvsPCgen_OutputFcn, ...
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
function varargout = PCvsPCgen_OutputFcn(~, ~, handles)
varargout{1} = handles.output;
end

% --- Executes JUST BEFORE OPENING
function PCvsPCgen_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.

% Choose default command line output for PCvsPCgen
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

set(handles.figure1,'toolbar','figure');
set(handles.figure1,'menubar','figure');

passedInHandles = varargin{1};
handles.stack = passedInHandles.stack;
handles.wnn = passedInHandles.wnn;
handles.rawDatasets = passedInHandles.rawDatasetsPassed;
handles.colorMatrix = passedInHandles.colorMatrix;
handles.colorMatrix2 = passedInHandles.colorMatrix2;
handles.colorMatrixRGBtrip = passedInHandles.colorMatrixRGBtrip;
handles.rootPath = passedInHandles.rootPath;
handles.standardsPath = passedInHandles.standardsPath;
handles.coeff = passedInHandles.coeff;
scores = passedInHandles.scores;

handles.scores = passedInHandles.scores;
handles.pcvars = passedInHandles.pcvars;

handles.cropAbove = passedInHandles.cropAbove;
handles.cropBelow = passedInHandles.cropBelow;

guidata(hObject,handles); % and then save them

%% set up the PCA dropdown boxes with the possible linkage distance metrics
% and methods
% default: 'ward' and 'euclidean'
methods = {'ward','average','centroid','complete','median','single','weighted'};
metrics = {'euclidean','seuclidean','cityblock','minkowski','chebychev','mahalanobis','cosine','correlation','spearman','hamming','jaccard'};
set(handles.linkageMetricsMenu,'String',metrics)
set(handles.linkageMethodsMenu,'String',methods)

% set the slider range and step size
numSteps = 10;
set(handles.maxClustersSlider, 'Min', 1);
set(handles.maxClustersSlider, 'Max', numSteps);
set(handles.maxClustersSlider, 'Value', 5);
set(handles.maxClustersSlider, 'SliderStep',[1/(numSteps-1),1/(numSteps-1)]);


%% plot the scores boxplot

axes(handles.scoresBoxplot)
cla; hold on;
set(gca,'Tag','scoresBoxplot')
xlabel('PC number')
ylabel('scores spread')

g = size(scores,2);
if g >= 10
    boxplot(scores(:,1:10),'orientation','vertical')
else
    boxplot(scores(:,1:g),'orientation','vertical')
end

%% plot the global Mean
axes(handles.globalMeanAxes)
cla; hold on;
set(gca,'Tag','globalMeanAxes')
shadedErrorBar(handles.wnn,mean(handles.stack,1),std(handles.stack,1),'k',2)
plot(handles.wnn,mean(handles.stack,1),'k','LineWidth',1)
set(gca,'xlim',[handles.cropBelow handles.cropAbove])
ylabel('Intensity (A.U.)');
xlabel('Rel. wavenumber (cm^-^1)')

%% clear the cluster graphs if necessary (e.g. window was already open)
% clear the cluster figures and tables
for z=1:3
    field = ['clusterAxes',num2str(z)];
    axes(handles.(field))
    cla
    set(gca,'Tag',field)
    
    field = ['clusterMenu',num2str(z)];
    set(handles.(field),'String','')
    set(handles.(field),'Value',[])
end


%% set up the drop box boxes that list each PC
numOfScores = size(scores,1)-1;
PCcatGen = @(number) ['PC',number,''];

listOfScores = cell(numOfScores,1);
for x=1:numOfScores
    listOfScores{x} = PCcatGen(num2str(x));
end

set(handles.PC1popupmenu,'String',listOfScores);
set(handles.PC1popupmenu,'Value',1);
set(handles.fromPCMenu,'String',listOfScores);
set(handles.fromPCMenu,'Value',1);
set(handles.toPCMenu,'String',listOfScores);
set(handles.toPCMenu,'Value',10);
set(handles.PC2popupmenu,'String',listOfScores);
set(handles.PC2popupmenu,'Value',2);

%% Plot the average spectra and first 4 PC coeffs
g = size(handles.scores,2);
for n=1:5
    field = ['leftPC',num2str(n)];
    axes(handles.(field))
    cla
    hold on
    plotCoeffs(n,field);
    hline(0)
    set(gca,'xlim',[handles.cropBelow handles.cropAbove])
    
end

for n=1:4
    field = ['bottomPC',num2str(n)];
    axes(handles.(field))
    cla
    hold on
    plotCoeffs(n,field);
    hline(0)
    set(gca,'xlim',[handles.cropBelow handles.cropAbove])
end

    function plotCoeffs(x,nameSent)
        if x > g
            return
        end
        plot(handles.wnn,handles.coeff(:,x),handles.colorMatrix{x})
        set(gca, 'YTickLabel',[]);
        if x ~= 6
            set(gca, 'XTickLabel',[])
        end
        set(gca,'Tag',nameSent)
    end

%% Plot PC scores
colergen1 = @(color1,color2,color3,text) ['<html><table border=0 width=400 font color=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];
colergen2 = @(color1,color2,color3,text) ['<html><table border=0 width=400 bgcolor=rgb(',color1,',',color2,',',color3,',',')><TR><TD><b>',text,'</b></TD></TR> </table></html>'];

colorsss = cellfun(@(x) x*255,handles.colorMatrix2','un',0); % convert to rgb style for html

for l=1:size(handles.rawDatasets,1)
    data{l,1} = colergen2(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),'');
    data{l,2} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasets{l,2}); % set the #samples column
    data{l,3} = colergen1(num2str(colorsss{l}(1)),num2str(colorsss{l}(2)),num2str(colorsss{l}(3)),handles.rawDatasets{l,1}); % set the dataset column
end

set(handles.samplesTable, 'Data',data)

dummyYvaluesOnes = ones(size(scores,1), 1);
dummyYvalues = dummyYvaluesOnes;
textSpacing = 0.1;

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

sampleNumbers = num2str(sampleNumbers);
numOfSamples = sum(numOfSamples); % then sum them up

handles.sampleNumbers = sampleNumbers;
handles.numOfSamples = numOfSamples;
handles.colorMatrixFromData = colorMatrixFromData;
handles.dummyYvaluesText = dummyYvaluesText;
handles.dummyYvalues = dummyYvalues;
guidata(hObject, handles);

axes(handles.mainPCplot)
cla
set(handles.mainPCplot,'Tag','mainPCplot')
redrawPlot(handles,2,1)

guidata(hObject, handles);

end

%% --- Re draw plots based on button or dropdown box selection
function PC1popupmenu_Callback(~, ~, handles)

currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

redrawPlot(handles,x,y)
end
function PC2popupmenu_Callback(~, ~, handles)
currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

redrawPlot(handles,x,y)
end
function leftNextbutton_Callback(~, ~, handles)

currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2})+1;

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

redrawPlot(handles,x,y)
end
function bottomNextButton_Callback(~, ~, handles)

currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2})+1;

redrawPlot(handles,x,y)
end
function PC1vsPC2_Callback(~, ~, handles)
x=2;
y=1;
redrawPlot(handles,x,y)
end
function PC1vsPC3_Callback(~, ~, handles)
x=3;
y=1;
redrawPlot(handles,x,y)
end
function PC2vsPC3_Callback(~, ~, handles)
x=3;
y=2;
redrawPlot(handles,x,y)
end

function redrawPlot(handles,x,y)

h = gcf;
allLines = findobj(h,'type', 'line');
for a=1:length(allLines)
    allLines(a).LineWidth = 1;
end

allPlots = findobj(h,'type', 'axes');
for a=1:length(allPlots)
    allPlots(a).Color = rgb('White');
end

extendPos = 1.2;
extendNeg = 1/extendPos;

axes(handles.mainPCplot)
cla
rotate3d off
view(2)
set(handles.mainPCplot, 'Tag','mainPCplot')
scatter(handles.scores(:,x),handles.scores(:,y),30,handles.colorMatrixFromData, 'filled')

for z=1:length(handles.colorMatrixFromData)
    handles.plotText = text(handles.scores(z,x),handles.scores(z,y),handles.sampleNumbers(z,:),'Color',handles.colorMatrixFromData(z,:),'FontSize',12);
end

guidata(handles.mainPCplot, handles);

hold on
xlabelgen = strcat('PC',num2str(x));
ylabelgen = strcat('PC',num2str(y));
xlabel(xlabelgen,'FontSize',14,'FontWeight', 'bold','Color',handles.colorMatrixRGBtrip{x});
ylabel(ylabelgen,'FontSize',14,'FontWeight', 'bold','Color',handles.colorMatrixRGBtrip{y});
set(handles.PC1popupmenu, 'Value', y)
set(handles.PC2popupmenu, 'Value', x)

PClabelGenY = strcat('leftPC',num2str(y));
PClabelGenX = strcat('bottomPC',num2str(x));

h = findobj(0,'Tag', PClabelGenY);
try
    plot = findobj(h,'type','line');
    plot.LineWidth = 2;
    h.Color = rgb('LightGray');
catch
end

h2 = findobj(0,'Tag', PClabelGenX);
try
    plot = findobj(h2,'type','line');
    plot.LineWidth = 2;
    h2.Color = rgb('LightGray');
catch
end

% Reset limits of the graph a bit to contain all the data
Xmax = max(handles.scores(:,x)); % find the mins and maxs for the chosen data
Xmin = min(handles.scores(:,x));
Ymax = max(handles.scores(:,y));
Ymin = min(handles.scores(:,y));

% extend them out by extendPos (or
if Xmax > 0
    Xmax = extendPos*Xmax;
else
    Xmax = extendNeg*Xmax;
end

if Ymax > 0
    Ymax = extendPos*Ymax;
else
    Ymax = extendNeg*Ymax;
end

if Xmin > 0
    Xmin = extendNeg*Xmin;
else
    Xmin = extendPos*Xmin;
end

if Ymin > 0
    Ymin = extendNeg*Ymin;
else
    Ymin = extendPos*Ymin;
end

% save the new limits
XLim = [Xmin Xmax];
YLim = [Ymin Ymax];

vline(0)
hline(0)
set(gca,'ButtonDownFcn',@(hObject,eventdata)PCvsPCgen('mainPCplot_ButtonDownFcn',hObject,eventdata,guidata(hObject)))

% and apply them to the main 2D scatter plot
set(handles.mainPCplot,'XLim',XLim);
set(handles.mainPCplot,'YLim',YLim);

PCname = ['fig-PC',num2str(y),'vsPC',num2str(x)];
set(handles.datasetPCNameEdit,'String',PCname)

end

%% LINKAGE and CLUSTER
function applyLinkagebutton_Callback(~, ~, handles)

currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

plotLinkage(handles,x,y,0.2)

set(handles.overlayCheckbox,'Enable','on')
set(handles.plotDendrogramButton,'Enable','on')
set(handles.fitClustersButton,'Enable','on')

end

function plotLinkage(handles,w,z,alphaValue)

%% LINKAGE AND CLUSTER ANALYSIS
% grab the 'from' and 'to' range for PC linkage user input dropdown boxes
f = cellstr(get(handles.fromPCMenu,'String'));
f = strsplit(f{get(handles.fromPCMenu,'Value')},'C');
f = str2double(f{2});

g = cellstr(get(handles.toPCMenu,'String'));
g = strsplit(g{get(handles.toPCMenu,'Value')},'C');
g = str2double(g{2});

% also grab the user input metric and method
contents = cellstr(get(handles.linkageMethodsMenu,'String'));
linkageMethod = contents{get(handles.linkageMethodsMenu,'Value')};
contents2 = cellstr(get(handles.linkageMetricsMenu,'String'));
linkageMetric = contents2{get(handles.linkageMetricsMenu,'Value')};

% LINKAGE
linKage = linkage(handles.scores(:,f:g),linkageMethod,linkageMetric);
handles.linKage = linKage; % and save it

% use the pdist in the desired PC range to optimally re-order the leafs in
% the dendrogram
D = pdist(handles.scores(:,f:g));
leafOrder = optimalleaforder(linKage,D);
handles.leafOrder = leafOrder; % and save it

% CLUSTER
% grab the maximum number of clusters for the
maxClustN = str2double(get(handles.maxClusterText,'String'));
cluSter = cluster(linKage,'maxclust',maxClustN);
handles.cluSter = cluSter; % and save it

guidata(gcbo, handles); % didn't pass hObject, so use gcbo here instead to save handles

%% GRAPH CLUSTER LINKAGE

axes(handles.mainPCplot)
cla
set(handles.mainPCplot,'Tag','mainPCplot')
redrawPlot(handles,w,z)


% this 'for' loop
% 1. looks through the cluSter matrix and pulls out the % coordinates for
% each point on the current PC vs PC axes
% 2.
for x=1:max(cluSter) % for each cluster
    % 1. grab the row numbers for the current cluster
    indices = find(cluSter == x);
    % 2. then store the vertices of each data point in this cluster
    vertices = zeros(length(indices),2); % initialize for loop
    for y=1:length(indices)
        stax(x).stack(y,:) = handles.stack(indices(y),:);
        scores1 = handles.scores(:,z);
        scores2 = handles.scores(:,w);
        vertices(y,1) = scores1(indices(y));
        vertices(y,2) = scores2(indices(y));
    end
    
    % 3. draw a filled polygon (or empty if user input boundaryOnlyCheckbox == 1)
    a = vertices(:,2);
    b = vertices(:,1);
    k = boundary(a,b); % draw a boundary around all the vertices in this cluster
    hold on;
    
    if length(a)>=5
        ellipse_t = fit_ellipse(a,b,handles.mainPCplot,handles.colorMatrix{x},alphaValue);
    else
        ellipse_t = [];
    end
    if isempty(ellipse_t)
        plot(a(k),b(k),handles.colorMatrix{x},'LineWidth',1,'LineStyle','-.');
    elseif strcmp(ellipse_t.status,'Hyperbola found')
        plot(a(k),b(k),handles.colorMatrix{x},'LineWidth',1,'LineStyle','-.');
    end
    
    tf = get(handles.boundaryOnlyCheckbox,'Value'); % return bool for checkbox state
    if tf == 1
        plot(a(k),b(k),handles.colorMatrix{x},'LineWidth',1); % plot the empty boundary
    end
end

handles.stax = stax;
guidata(gcbo, handles);

% after initializing the linkage for the first time, the info can be used
% to draw a dendrogram
set(handles.plotDendrogramButton,'Enable','on')
set(handles.clearClustersFromGraphButton,'Enable','on')

% graph the average spectra difference from the mean for each cluster
% defaulted to plot the first three in the ClustersPanel

% This chunk of code generates the clusterMenu tables for the ClustersPanel
% and then graphs the first three clusters
colergen1 = @(color1,color2,color3,text) ['<html><font color=rgb(',color1,',',color2,',',color3,',',')><b>',text,'</b></font></html>'];
colorsss = cellfun(@(x) x*255,handles.colorMatrix2','un',0); % convert to rgb style for html
clusterNames = cell(max(cluSter),1);
for x=1:max(cluSter)
    thisClusterName = ['Cluster ',num2str(x)];
    clusterNames(x) = cellstr(colergen1(num2str(colorsss{x}(1)),num2str(colorsss{x}(2)),num2str(colorsss{x}(3)),thisClusterName));
end

handles.clusterNames = clusterNames;
guidata(gcbo, handles); % save to pass to fitSpectra.m

% if there is less than 3 clusters specified (or similarly less than 2),
% the 2nd/3rd axes and listboxes are cleared
set(handles.clusterMenu1,'String',clusterNames) % set up clusterMenu1
set(handles.clusterMenu1,'Value',1)
if max(cluSter) == 2
    set(handles.clusterMenu2,'String',clusterNames) % set up clusterMenu1
    set(handles.clusterMenu2,'Value',2)
    set(handles.clusterMenu3,'String','') % set up clusterMenu1
    set(handles.clusterMenu3,'Value',[])
elseif max(cluSter) >= 3
    set(handles.clusterMenu2,'String',clusterNames) % set up clusterMenu1
    set(handles.clusterMenu2,'Value',2)
    set(handles.clusterMenu3,'String',clusterNames) % set up clusterMenu1
    set(handles.clusterMenu3,'Value',3)
elseif max(cluSter) == 1
    set(handles.clusterMenu2,'String','') % set up clusterMenu1
    set(handles.clusterMenu2,'Value',[])
    set(handles.clusterMenu3,'String','') % set up clusterMenu1
    set(handles.clusterMenu3,'Value',[])
end

if max(cluSter) <= 3
    numClusters = max(cluSter);
else
    numClusters = 3;
end

handles.numClusters = numClusters;
guidata(gcbo, handles); % save to pass to fitSpectra.m

axes(handles.clusterAxes1); cla; % clear the figures
axes(handles.clusterAxes2); cla;
axes(handles.clusterAxes3); cla;

for x=1:numClusters
    field = ['clusterAxes',num2str(x)];
    axes(handles.(field))
    set(gca,'Tag',field)
    hold on
    hline(0)
    tf = get(handles.overlayCheckbox,'Value');
    if tf == 1
        shadedErrorBar(handles.wnn,mean(handles.stack,1),std(handles.stack,1),'k',2)
        plot(handles.wnn,mean(handles.stack,1),'k','LineWidth',1,'LineStyle','--')
    end
    plot(handles.wnn,(mean(stax(x).stack,1)-mean(handles.stack,1)),handles.colorMatrix{x},'LineWidth',2)
    set(gca,'xlim',[handles.cropBelow handles.cropAbove])
end



end

function plotDendrogramButton_Callback(~, ~, handles)

try
    figure(4)
    hold off
    H = dendrogram(handles.linKage, 0,'Reorder',handles.leafOrder,'orientation','left','ColorThreshold','default');
catch
    messageBox('Must initialize linkage before plotting', 'Error','error')
end


end

function maxClusterText_Callback(hObject, ~, ~)
% this callback makes sure no value between n and l can be entered
n = 1; % change these limits if desired
l = 10;
maxClust = str2double(get(hObject,'String')); %returns contents of maxClusterText as a double

if maxClust > l
    set(hObject,'String','10')
elseif maxClust < n
    set(hObject,'String','1')
else
    set(hObject,'String',num2str(maxClust))
end

end
function maxClustersSlider_Callback(hObject, ~, handles)

% this rounds to the nearest integer
val=round(hObject.Value);
hObject.Value=val;
set(handles.maxClusterText,'String',num2str(val))

end

function clearClustersFromGraphButton_Callback(~, ~, handles)
currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

redrawPlot(handles,x,y)

% clear the cluster figures and tables
for z=1:3
    field = ['clusterAxes',num2str(z)];
    axes(handles.(field))
    cla
    
    set(gca,'Tag',field)
    
    field = ['clusterMenu',num2str(z)];
    set(handles.(field),'String','')
    set(handles.(field),'Value',[])
end

end

% --- Executes on button press in exportPlotButton.
function exportPlotButton_Callback(~, ~, handles)
% exports the mainPlotPC to a vector pdf for further manipulation
% plotLinkage(handles,x,y,1)
% 1. grab the user-input or default name from the datasetPCNameEdit box
fileName = [get(handles.datasetPCNameEdit,'String'),'.pdf']; % and append .pdf

% '-q101' - sets the compression to lossless
%'-cmyk' - sets the appropriate color scheme
% '-painters' - helps with dotted lines, etc.
export_fig(handles.mainPCplot, fileName,'-q101','-cmyk','-painters');

%  plotLinkage(handles,x,y,0.2)
fprintf([fileName,' exported to current directory...\n']);
end


%% hanging CreateFcn's
% most of these just set up the white bckgnd for listboxes and popupmenus
function listbox1_CreateFcn(hObject, ~, ~) %#ok<*DEFNU>
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function linkageMetricsMenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function linkageMethodsMenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function maxClusterText_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function maxClustersSlider_CreateFcn(hObject, ~, ~)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
end
function fromPCMenu_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function toPCMenu_CreateFcn(hObject, ~, ~)
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
function clusterMenu1_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function clusterMenu2_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function clusterMenu3_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
function datasetPCNameEdit_CreateFcn(hObject, ~, ~)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

%% Clusters Panel

% --- Executes on selection change in clusterMenu1.
function clusterMenu1_Callback(hObject, ~, handles)
try
    contents = cellstr(get(hObject,'String')); %returns testy contents as cell array
    selCluster = strsplit(contents{get(hObject,'Value')},{'Cluster ','</b>'}); %returns selected item from testy
    numCluster = str2double(selCluster{2});
    plotNewCluster(1,numCluster,handles)
catch
    fprintf('Calculate linkage before analyzing clusters\n');
end

end

% --- Executes on selection change in clusterMenu2.
function clusterMenu2_Callback(hObject, ~, handles)
try
    contents = cellstr(get(hObject,'String')); %returns testy contents as cell array
    selCluster = strsplit(contents{get(hObject,'Value')},{'Cluster ','</b>'}); %returns selected item from testy
    numCluster = str2double(selCluster{2});
    plotNewCluster(2,numCluster,handles)
catch
    fprintf('Calculate linkage before analyzing clusters\n');
end

end

% --- Executes on selection change in clusterMenu3.
function clusterMenu3_Callback(hObject, ~, handles)
try
    contents = cellstr(get(hObject,'String')); %returns testy contents as cell array
    selCluster = strsplit(contents{get(hObject,'Value')},{'Cluster ','</b>'}); %returns selected item from testy
    numCluster = str2double(selCluster{2});
    plotNewCluster(3,numCluster,handles)
catch
    fprintf('Calculate linkage before analyzing clusters\n');
end

end

function plotNewCluster(fig,n,handles)

field = ['clusterAxes',num2str(fig)];
axes(handles.(field))
cla
set(gca,'Tag',field)
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
% --- Executes on button press in overlayCheckbox.
function overlayCheckbox_Callback(~, ~, handles)
currentListOfScores = cellstr(get(handles.PC1popupmenu,'String'));
currentSelection = get(handles.PC1popupmenu,'Value');
split = strsplit(currentListOfScores{currentSelection},'C');
y = str2double(split{2});

currentListOfScores2 = cellstr(get(handles.PC2popupmenu,'String'));
currentSelection2 = get(handles.PC2popupmenu,'Value');
split2 = strsplit(currentListOfScores2{currentSelection2},'C');
x = str2double(split2{2});

plotLinkage(handles,x,y,0.2)
end

function datasetPCNameEdit_Callback(~, ~, ~)
end

%% To implement...

% --- Executes on mouse press over axes background.
function mainPCplot_ButtonDownFcn(~, ~, ~)
% want to develop code here that is capable of selecting points for
% deletion
a = get(gca,'currentpoint');
end

% --- Executes on button press in fitClustersButton.
function fitClustersButton_Callback(~, ~, handles)
fitSpectra(handles)
end


% --- Executes when selected cell(s) is changed in samplesTable.
function samplesTable_CellSelectionCallback(hObject, eventdata, handles)
% upon selection of a dataset in the listbox, change
% 1. text to bold and
% 2. marker type to x's

% first, clear any previous text changes

plotText = flip(findobj(gcf,'type', 'text')); % returns backwards so need to flip
for x=1:length(plotText)
    plotText(x).FontWeight = 'normal';
    plotText(x).FontSize = 12;
end

tableData = get(hObject, 'data');

% use eventdata to grab the appropriate row
selRow = eventdata.Indices(1);

% get the range of the sample numbers for the selected dataset
samplesBeforeOurRow = 0;
if selRow ~= 1
    for x=1:(selRow-1)
        split1 = strsplit(tableData{x,2},'<b>');
        split2 = strsplit(split1{2},'</b>');
        numSamples = str2double(split2{1});
        samplesBeforeOurRow = samplesBeforeOurRow + numSamples;
    end
end

% get numSamples in selected dataset
split1 = strsplit(tableData{selRow,2},'<b>');
split2 = strsplit(split1{2},'</b>');
numSamples = str2double(split2{1});

sampleStart = samplesBeforeOurRow + 1;
sampleEnd = samplesBeforeOurRow + numSamples;

plotText = flip(findobj(gcf,'type', 'text')); % returns the text backwards so need to flip

for x=sampleStart:sampleEnd
    plotText(x).FontWeight = 'bold';
    plotText(x).FontSize = 20;
end

end
