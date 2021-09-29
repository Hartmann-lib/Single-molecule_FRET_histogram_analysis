function varargout = FitXGauss(varargin)
% Graphical user interface to fit multiple Gaussian distributions to 
% a cumulated FRET efficiency histogram collected from individual (weighted) 
% TIRF-FRET traces

% 2017-20-11 - Andreas Hartmann 

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FitXGauss_OpeningFcn, ...
                   'gui_OutputFcn',  @FitXGauss_OutputFcn, ...
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

function FitXGauss_OpeningFcn(hObject, eventdata, handles, varargin)

global ButtonLayout;
global ButtonLayout2;
global ButtonLayout3;
global barColor;

addpath(genpath('scripts'));

handles.output = hObject;

guidata(hObject, handles);
set(handles.pushbutton_FitGauss,'Enable','off');
set(handles.pushbuttonExport,'Enable','off');
set(handles.pushbuttonInput,'Enable','off');
set(handles.BootSample,'Enable','off');
set(handles.BootFits,'Enable','off');
set(handles.Probs,'Enable','off');
set(handles.checkWeight,'Enable','off');
set(handles.checkNorm,'Enable','off');
set(handles.popupGauss,'Enable','off');
set(handles.checkSumFit,'Enable','off');
set(handles.pushbutton_FitGauss,'CData',[]);
set(handles.axesFRET,'Box','on');
xlim(handles.axesFRET,[-0.1 1.1]);
xlabel(handles.axesFRET,'FRET efficiency');
ylabel(handles.axesFRET,'Counts');

ButtonLayout=importdata('Gauss_small.tif');
ButtonLayout2=importdata('download2.tif');
ButtonLayout3=importdata('target.tif');

set(handles.uitableGauss,'Data',[]);

barColor=[0.75 0.75 0.75];

function varargout = FitXGauss_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.output;

function checkWeight_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global weightFRET;
    global sumFRET;
    global barColor;
    
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    set(handles.BootSample,'Enable','off');
    set(handles.BootFits,'Enable','off');
    set(handles.Probs,'Enable','off');
    set(handles.pushbuttonExport,'CData',[]);    
    set(handles.pushbuttonExport,'Enable','off');
    set(handles.uitableResults,'Data',[]);
    
    cla(handles.axesFRET);
    
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    
    if boolWeights
        
        showFRET=weightFRET./(1-boolNorm+boolNorm*sum(weightFRET));
    else
        
        showFRET=sumFRET./(1-boolNorm+boolNorm*sum(sumFRET));
    end
        
    b1=bar(handles.axesFRET,edgesFRET,showFRET,'histc');
    set(b1,'FaceColor',barColor);
    xlim(handles.axesFRET,[-0.1 1.1]);
    xlabel(handles.axesFRET,'FRET efficiency');
    
    if boolWeights
               
        if boolNorm
            ylabel(handles.axesFRET,'\itp\rm(\Delta\itE\rm)');
        else
            ylabel(handles.axesFRET,'Number of molecules');
        end
    else
        if boolNorm
            ylabel(handles.axesFRET,'\itp\rm(\Delta\itE\rm)');
        else
            ylabel(handles.axesFRET,'Number of frames');
        end
    end
    
    popupGauss_Callback(hObject, eventdata, handles);
    axis 'auto y';

function loadHist_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global weightFRET;
    global sumFRET;
    global barColor;
    global molFRET;
    
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    
    cla(handles.axesFRET);    
    
    set(handles.popupGauss,'Value',1);
    set(handles.popupGauss,'Enable','off');
    popupGauss_Callback(hObject, eventdata, handles);
    set(handles.checkWeight,'Enable','off');
    set(handles.checkNorm,'Enable','off');
    set(handles.pushbuttonExport,'Enable','off');
    set(handles.pushbuttonExport,'CData',[]);
    
    boolWeights=boolean(get(handles.checkWeight,'Value'));

    % Select data files
    [filenames,pathname]=uigetfile('*.hist2','MultiSelect','on', 'Histogram Data:');
    
    set(handles.figure1,'Name',['FitXGauss - ' pathname]);
    
    if ~isempty(filenames)
    
        if ischar(filenames)

            fileArr=cell(1,1);
            fileArr{1}=filenames;
        else

            fileArr=filenames;    
        end

        % Number of files
        numFiles=length(fileArr);

        h=waitbar(0,'Please wait...');

        molFRET={[]};
        
        for iter=1:numFiles

            currData=dlmread([pathname char(fileArr(iter))]);

            if iter==1

                edgesFRET=currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,1);
                sumFRET=currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2);
                weightFRET=currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2)./sum(currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2));
                molFRET{iter,1}=currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2);
                molFRET{iter,2}=sum(currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2));
            else

                sumFRET=sumFRET+currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2);
                weightFRET=weightFRET+currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2)./sum(currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2));
                molFRET{iter,1}=currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2);
                molFRET{iter,2}=sum(currData(currData(:,1)>=-0.1&currData(:,1)<=1.1,2));
            end

            waitbar(iter/numFiles);
        end

        close(h);

        if boolWeights
            
            showFRET=weightFRET./(1-boolNorm+boolNorm*sum(weightFRET));
        else
            
            showFRET=sumFRET./(1-boolNorm+boolNorm*sum(sumFRET));
        end

        b1=bar(handles.axesFRET,edgesFRET,showFRET,'histc');
        set(b1,'FaceColor',barColor);
        xlim(handles.axesFRET,[-0.1 1.1]);
        xlabel(handles.axesFRET,'FRET efficiency');
        
        if boolWeights
            
            if boolNorm
                ylabel(handles.axesFRET,'\itp\rm(\Delta\itE\rm)');
            else
                ylabel(handles.axesFRET,'Number of molecules');
            end
        else
            if boolNorm
                ylabel(handles.axesFRET,'\itp\rm(\Delta\itE\rm)');
            else
                ylabel(handles.axesFRET,'Number of frames');
            end
        end

        set(handles.numMolecules,'String',[num2str(numFiles) ' Molecules']);

        set(handles.popupGauss,'Enable','on');
        set(handles.checkWeight,'Enable','on');
        set(handles.checkNorm,'Enable','on');
        
        set(handles.editSamples,'String',num2str(numFiles));
        set(handles.editIterations,'String',num2str(numFiles));
    end

function File_Callback(hObject, eventdata, handles)

function popupGauss_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global weightFRET;
    global sumFRET;    
    global ButtonLayout;
    global ButtonLayout3;

    numGauss=get(handles.popupGauss,'Value')-1;
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    
    if boolWeights
               
        showFRET=weightFRET./(1-boolNorm+boolNorm*sum(weightFRET));  
    else
        
        showFRET=sumFRET./(1-boolNorm+boolNorm*sum(sumFRET));
    end
    
    maxFRET=max(showFRET);
    
    if numGauss>0
        
        switch numGauss
            case 1
                uiData=[round(maxFRET) 0 max(showFRET) sum(edgesFRET.*showFRET)/sum(showFRET) 0 1 0.1 0 1];
            case 2
                uiData=[round(maxFRET/numGauss) 0 maxFRET 1/(1+numGauss) 0 1/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 2/(1+numGauss) 1/numGauss 1 0.1 0 1];                
            case 3
                uiData=[round(maxFRET/numGauss) 0 maxFRET 1/(1+numGauss) 0 1/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 2/(1+numGauss) 1/numGauss 2/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 3/(1+numGauss) 2/numGauss 3/numGauss 0.1 0 1];                
            case 4
                uiData=[round(maxFRET/numGauss) 0 maxFRET 1/(1+numGauss) 0 1/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 2/(1+numGauss) 1/numGauss 2/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 3/(1+numGauss) 2/numGauss 3/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 4/(1+numGauss) 3/numGauss 4/numGauss 0.1 0 1];                
            otherwise
                uiData=[round(maxFRET/numGauss) 0 maxFRET 1/(1+numGauss) 0 1/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 2/(1+numGauss) 1/numGauss 2/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 3/(1+numGauss) 2/numGauss 3/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 4/(1+numGauss) 3/numGauss 4/numGauss 0.1 0 1;round(maxFRET/numGauss) 0 maxFRET 5/(1+numGauss) 4/numGauss 5/numGauss 0.1 0 1];                
        end

        set(handles.uitableGauss,'Data',uiData);
        set(handles.pushbutton_FitGauss,'Enable','on');
        set(handles.pushbutton_FitGauss,'CData',ButtonLayout);
        set(handles.pushbuttonInput,'Enable','on');
        set(handles.pushbuttonInput,'CData',ButtonLayout3);
    else
        set(handles.uitableGauss,'Data',[]);
        set(handles.pushbutton_FitGauss,'Enable','off');
        set(handles.pushbutton_FitGauss,'CData',[]);
        set(handles.pushbuttonInput,'Enable','off');
        set(handles.pushbuttonInput,'CData',[]);
    end
    
function popupGauss_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function pushbutton_FitGauss_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global weightFRET;
    global sumFRET;
    global barColor;
    global ButtonLayout2;
    global edgesFine;
    global cfunX;
    global cfunXMin;
    global cfunXMax;
    global molFRET;
    global hTotal;
    global bootFRET;
    global bootFit;
    global bootFitAll;
    global results2EXP;
    
    % load user inputs
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    boolTotal=boolean(get(handles.checkSumFit,'Value'));
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    numSamples=str2double(get(handles.editSamples,'String'));
    numIterations=str2double(get(handles.editIterations,'String'));
    numFiles=size(molFRET,1);
    
    % load colormap for fit curves
    colorLines=colormap('lines');
    
    % draw molecuels with replacement
    if boolWeights
        
        randIndex=randi(numFiles,numSamples,numIterations);
    else
        
        randIndex=randi_hist((1:numFiles)',cell2mat(molFRET(:,2)),numSamples,numIterations);
    end
    
    bootFRET=zeros(length(edgesFRET),numIterations); 
    
    for iter1=1:numIterations
        
        for iter2=1:numSamples

            currFRET=molFRET{randIndex(iter2,iter1),1};

            if boolWeights

                bootFRET(:,iter1)=bootFRET(:,iter1)+currFRET./sum(currFRET);
            else
                bootFRET(:,iter1)=bootFRET(:,iter1)+currFRET;
            end
        end
        
        bootFRET(:,iter1)=bootFRET(:,iter1)./numSamples.*numFiles;
        
        if boolNorm
            
            bootFRET(:,iter1)=bootFRET(:,iter1)./sum(bootFRET(:,iter1));
        end
    end
    
    % load fit settings
    uiData = get(handles.uitableGauss,'Data');
   
    startD=[];
    upperD=[];
    lowerD=[];
    for iter=1:size(uiData,1)
        
        startD=[startD uiData(iter,1) uiData(iter,4) uiData(iter,7)];
        lowerD=[lowerD uiData(iter,2) uiData(iter,5) uiData(iter,8)];
        upperD=[upperD uiData(iter,3) uiData(iter,6) uiData(iter,9)];
    end
      
    % Fit measured histogram
    if boolWeights
               
        showFRET=weightFRET./(1-boolNorm+boolNorm*sum(weightFRET));  
    else
        
        showFRET=sumFRET./(1-boolNorm+boolNorm*sum(sumFRET));
    end    
    
    % BOOTSTRAPPING
    %%%%%%%%%%%%%%%
    paramsMEAS=zeros(4*size(uiData,1),1);
    
    [cfunMEAS,gof,output]=fit(edgesFRET,showFRET,['gauss' num2str(size(uiData,1))],'Start',startD,'Lower',lowerD,'Upper',upperD);
        
    % Plot fit result
    cla(handles.axesFRET);
    hold(handles.axesFRET,'on');
    b1=bar(handles.axesFRET,edgesFRET,showFRET,'histc');
    set(b1,'FaceColor',barColor);
    xlim(handles.axesFRET,[-0.1 1.1]);
    xlabel(handles.axesFRET,'FRET efficiency');
    
    edgesFine=(-0.1:0.001:1.1);
    cfunX{1}=cfunMEAS.a1.*exp(-((edgesFine-cfunMEAS.b1)./cfunMEAS.c1).^2);
    
    % html generation for colored rows
    colergen = @(color,text) ['<html><table border=0 width=400 color=#FFFFFF bgcolor=',color,'><TR><TD>',text,'</TD></TR> </table></html>'];
      
    % Array of fit results
    sumArea=trapz(edgesFine,cfunMEAS(edgesFine));
    A1=cfunMEAS.a1*cfunMEAS.c1*sqrt(pi);
    
    p1{1}=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{1},'-k','LineWidth',2);
    set(p1{1},'Color',colorLines(1,:));
    
    if size(uiData,1)>1
        cfunX{2}=cfunMEAS.a2.*exp(-((edgesFine-cfunMEAS.b2)./cfunMEAS.c2).^2);
        p1{2}=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{2},'-k','LineWidth',2);
        set(p1{2},'Color',colorLines(2,:));
    end
    
    if size(uiData,1)>2
        cfunX{3}=cfunMEAS.a3.*exp(-((edgesFine-cfunMEAS.b3)./cfunMEAS.c3).^2);
        p1{3}=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{3},'-k','LineWidth',2);
        set(p1{3},'Color',colorLines(3,:));
   end
    
    if size(uiData,1)>3
        cfunX{4}=cfunMEAS.a4.*exp(-((edgesFine-cfunMEAS.b4)./cfunMEAS.c4).^2);
        p1{4}=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{4},'-k','LineWidth',2);
        set(p1{4},'Color',colorLines(4,:));
    end
    
    if size(uiData,1)>4
        cfunX{5}=cfunMEAS.a5.*exp(-((edgesFine-cfunMEAS.b5)./cfunMEAS.c5).^2);
        p1{5}=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{5},'-k','LineWidth',2);
        set(p1{5},'Color',colorLines(5,:));
    end

    hTotal=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunMEAS(edgesFine),'--k','LineWidth',2);
    
    if boolTotal
        set(hTotal,'Visible','on');
    else
        set(hTotal,'Visible','off');
    end
    
    % BOOTSTRAPPING
    %%%%%%%%%%%%%%%
    % Run gaussian fitting M times
    paramsBOOT=zeros(4*size(uiData,1),numIterations);
    
    h=waitbar(0,'Bootstrapping is running...');

    bootFit=zeros(size(uiData,1),length(edgesFine),numIterations);
    bootFitAll=zeros(length(edgesFine),numIterations);
        
    for iter=1:numIterations
                
       [cfun,gof,output]=fit(edgesFRET,bootFRET(:,iter),['gauss' num2str(size(uiData,1))],'Start',startD,'Lower',lowerD,'Upper',upperD);
       
%        figure();
%        bar(edgesFRET,bootFRET(:,iter),'hist');
%        hold on;
%        plot(cfun);
%        xlim([-0.1 1.1]);             
       
       bootFitAll(:,iter)=cfun(edgesFine);
       
       if size(uiData,1)==1
           
           bootFit(1,:,iter)=cfun.a1.*exp(-((edgesFine-cfun.b1)./cfun.c1).^2);          
           paramsBOOT(:,iter)=[cfun.a1 cfun.b1 cfun.c1 1];
       
       elseif size(uiData,1)==2
           
           bootFit(1,:,iter)=cfun.a1.*exp(-((edgesFine-cfun.b1)./cfun.c1).^2);
           bootFit(2,:,iter)=cfun.a2.*exp(-((edgesFine-cfun.b2)./cfun.c2).^2);
           A1=cfun.a1*cfun.c1*sqrt(pi);
           A2=cfun.a2*cfun.c2*sqrt(pi);
           paramsBOOT(:,iter)=[cfun.a1 cfun.b1 cfun.c1 A1/(A1+A2) cfun.a2 cfun.b2 cfun.c2 A2/(A1+A2)];
           
       elseif size(uiData,1)==3
           
           bootFit(1,:,iter)=cfun.a1.*exp(-((edgesFine-cfun.b1)./cfun.c1).^2);
           bootFit(2,:,iter)=cfun.a2.*exp(-((edgesFine-cfun.b2)./cfun.c2).^2);
           bootFit(3,:,iter)=cfun.a3.*exp(-((edgesFine-cfun.b3)./cfun.c3).^2);
           A1=cfun.a1*cfun.c1*sqrt(pi);
           A2=cfun.a2*cfun.c2*sqrt(pi);
           A3=cfun.a3*cfun.c3*sqrt(pi);
           paramsBOOT(:,iter)=[cfun.a1 cfun.b1 cfun.c1 A1/(A1+A2+A3) cfun.a2 cfun.b2 cfun.c2 A2/(A1+A2+A3) cfun.a3 cfun.b3 cfun.c3 A3/(A1+A2+A3)];

       elseif size(uiData,1)==4

           bootFit(1,:,iter)=cfun.a1.*exp(-((edgesFine-cfun.b1)./cfun.c1).^2);
           bootFit(2,:,iter)=cfun.a2.*exp(-((edgesFine-cfun.b2)./cfun.c2).^2);
           bootFit(3,:,iter)=cfun.a3.*exp(-((edgesFine-cfun.b3)./cfun.c3).^2);
           bootFit(4,:,iter)=cfun.a4.*exp(-((edgesFine-cfun.b4)./cfun.c4).^2);
           A1=cfun.a1*cfun.c1*sqrt(pi);
           A2=cfun.a2*cfun.c2*sqrt(pi);
           A3=cfun.a3*cfun.c3*sqrt(pi);
           A4=cfun.a4*cfun.c4*sqrt(pi);
           paramsBOOT(:,iter)=[cfun.a1 cfun.b1 cfun.c1 A1/(A1+A2+A3+A4) cfun.a2 cfun.b2 cfun.c2 A2/(A1+A2+A3+A4) cfun.a3 cfun.b3 cfun.c3 A3/(A1+A2+A3+A4) cfun.a4 cfun.b4 cfun.c4 A4/(A1+A2+A3+A4)];

       elseif size(uiData,1)==5
           
           bootFit(1,:,iter)=cfun.a1.*exp(-((edgesFine-cfun.b1)./cfun.c1).^2);
           bootFit(2,:,iter)=cfun.a2.*exp(-((edgesFine-cfun.b2)./cfun.c2).^2);
           bootFit(3,:,iter)=cfun.a3.*exp(-((edgesFine-cfun.b3)./cfun.c3).^2);
           bootFit(4,:,iter)=cfun.a4.*exp(-((edgesFine-cfun.b4)./cfun.c4).^2);
           bootFit(5,:,iter)=cfun.a5.*exp(-((edgesFine-cfun.b5)./cfun.c5).^2);
           A1=cfun.a1*cfun.c1*sqrt(pi);
           A2=cfun.a2*cfun.c2*sqrt(pi);
           A3=cfun.a3*cfun.c3*sqrt(pi);
           A4=cfun.a4*cfun.c4*sqrt(pi);
           A5=cfun.a5*cfun.c5*sqrt(pi);
           paramsBOOT(:,iter)=[cfun.a1 cfun.b1 cfun.c1 A1/(A1+A2+A3+A4+A5) cfun.a2 cfun.b2 cfun.c2 A2/(A1+A2+A3+A4+A5) cfun.a3 cfun.b3 cfun.c3 A3/(A1+A2+A3+A4+A5) cfun.a4 cfun.b4 cfun.c4 A4/(A1+A2+A3+A4+A5) cfun.a5 cfun.b5 cfun.c5 A5/(A1+A2+A3+A4+A5)];           
       end
       
       waitbar(iter/numIterations);
    end
    
    close(h);
    
    % delete old graphs
    for iter=1:length(p1)
        
        delete(p1{iter});
    end
    set(hTotal,'Visible','off');
    
    % Average fit params
    meanBOOT=mean(paramsBOOT,2);
    stdBOOT=std(paramsBOOT');

    lmrArray=lrSD(squeeze(bootFit(1,:,:)));
    cfunXMin{1}=lmrArray(:,2)-lmrArray(:,1);
    cfunX{1}=lmrArray(:,2);
    cfunXMax{1}=lmrArray(:,2)+lmrArray(:,3);
    
    actClr=sprintf('#%s%s%s',dec2hex(round(colorLines(1,1)*255)),dec2hex(round(colorLines(1,2)*255)),dec2hex(round(colorLines(1,3)*255)));
    
    results={colergen(actClr,num2str(meanBOOT(1))),colergen(actClr,num2str(stdBOOT(1))),colergen(actClr,num2str(meanBOOT(2))),colergen(actClr,num2str(stdBOOT(2))),colergen(actClr,num2str(meanBOOT(3)/sqrt(2))),colergen(actClr,num2str(stdBOOT(3)/sqrt(2))),colergen(actClr,num2str(meanBOOT(4))),colergen(actClr,num2str(stdBOOT(4)))};
    results2EXP=[meanBOOT(1),stdBOOT(1),meanBOOT(2),stdBOOT(2),meanBOOT(3)/sqrt(2),stdBOOT(3)/sqrt(2),meanBOOT(4),stdBOOT(4)];
     
    p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{1},'-k','LineWidth',2);
    set(p1,'Color',colorLines(1,:));
    p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMax{1},'--k','LineWidth',2);
    set(p1,'Color',colorLines(1,:));
    p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMin{1},'--k','LineWidth',2);
    set(p1,'Color',colorLines(1,:));
    
    if size(uiData,1)>1

        lmrArray=lrSD(squeeze(bootFit(2,:,:)));
        cfunXMin{2}=lmrArray(:,2)-lmrArray(:,1);
        cfunX{2}=lmrArray(:,2);
        cfunXMax{2}=lmrArray(:,2)+lmrArray(:,3);

        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{2},'-k','LineWidth',2);
        set(p1,'Color',colorLines(2,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMax{2},'--k','LineWidth',2);
        set(p1,'Color',colorLines(2,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMin{2},'--k','LineWidth',2);
        set(p1,'Color',colorLines(2,:));
        actClr=sprintf('#%s%s%s',dec2hex(round(colorLines(2,1)*255)),dec2hex(round(colorLines(2,2)*255)),dec2hex(round(colorLines(2,3)*255)));
        results(end+1,:)={colergen(actClr,num2str(meanBOOT(5))),colergen(actClr,num2str(stdBOOT(5))),colergen(actClr,num2str(meanBOOT(6))),colergen(actClr,num2str(stdBOOT(6))),colergen(actClr,num2str(meanBOOT(7)/sqrt(2))),colergen(actClr,num2str(stdBOOT(7)/sqrt(2))),colergen(actClr,num2str(meanBOOT(8))),colergen(actClr,num2str(stdBOOT(8)))};
        results2EXP=[results2EXP;meanBOOT(5),stdBOOT(5),meanBOOT(6),stdBOOT(6),meanBOOT(7)/sqrt(2),stdBOOT(7)/sqrt(2),meanBOOT(8),stdBOOT(8)];
    end
    
    if size(uiData,1)>2
        
        lmrArray=lrSD(squeeze(bootFit(3,:,:)));
        cfunXMin{3}=lmrArray(:,2)-lmrArray(:,1);
        cfunX{3}=lmrArray(:,2);
        cfunXMax{3}=lmrArray(:,2)+lmrArray(:,3);
        
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{3},'-k','LineWidth',2);
        set(p1,'Color',colorLines(3,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMin{3},'--k','LineWidth',2);
        set(p1,'Color',colorLines(3,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMax{3},'--k','LineWidth',2);
        set(p1,'Color',colorLines(3,:));
        actClr=sprintf('#%s%s%s',dec2hex(round(colorLines(3,1)*255)),dec2hex(round(colorLines(3,2)*255)),dec2hex(round(colorLines(3,3)*255)));
        results(end+1,:)={colergen(actClr,num2str(meanBOOT(9))),colergen(actClr,num2str(stdBOOT(9))),colergen(actClr,num2str(meanBOOT(10))),colergen(actClr,num2str(stdBOOT(10))),colergen(actClr,num2str(meanBOOT(11)/sqrt(2))),colergen(actClr,num2str(stdBOOT(11)/sqrt(2))),colergen(actClr,num2str(meanBOOT(12))),colergen(actClr,num2str(stdBOOT(12)))};
        results2EXP=[results2EXP;meanBOOT(9),stdBOOT(9),meanBOOT(10),stdBOOT(10),meanBOOT(11)/sqrt(2),stdBOOT(11)/sqrt(2),meanBOOT(12),stdBOOT(12)];
    end
    
    if size(uiData,1)>3
        
        lmrArray=lrSD(squeeze(bootFit(4,:,:)));
        cfunXMin{4}=lmrArray(:,2)-lmrArray(:,1);
        cfunX{4}=lmrArray(:,2);
        cfunXMax{4}=lmrArray(:,2)+lmrArray(:,3);
        
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{4},'-k','LineWidth',2);
        set(p1,'Color',colorLines(4,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMax{4},'--k','LineWidth',2);
        set(p1,'Color',colorLines(4,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMin{4},'--k','LineWidth',2);
        set(p1,'Color',colorLines(4,:));
        actClr=sprintf('#%s%s%s',dec2hex(round(colorLines(4,1)*255)),dec2hex(round(colorLines(4,2)*255)),dec2hex(round(colorLines(4,3)*255)));
        results(end+1,:)={colergen(actClr,num2str(meanBOOT(13))),colergen(actClr,num2str(stdBOOT(13))),colergen(actClr,num2str(meanBOOT(14))),colergen(actClr,num2str(stdBOOT(14))),colergen(actClr,num2str(meanBOOT(15)/sqrt(2))),colergen(actClr,num2str(stdBOOT(15)/sqrt(2))),colergen(actClr,num2str(meanBOOT(16))),colergen(actClr,num2str(stdBOOT(16)))};
        results2EXP=[results2EXP;meanBOOT(13),stdBOOT(13),meanBOOT(14),stdBOOT(14),meanBOOT(15)/sqrt(2),stdBOOT(15)/sqrt(2),meanBOOT(16),stdBOOT(16)];
    end
    
    if size(uiData,1)>4

        lmrArray=lrSD(squeeze(bootFit(5,:,:)));
        cfunXMin{5}=lmrArray(:,2)-lmrArray(:,1);
        cfunX{5}=lmrArray(:,2);
        cfunXMax{5}=lmrArray(:,2)+lmrArray(:,3);
        
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunX{5},'-k','LineWidth',2);
        set(p1,'Color',colorLines(5,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMax{5},'--k','LineWidth',2);
        set(p1,'Color',colorLines(5,:));
        p1=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),cfunXMin{5},'--k','LineWidth',2);
        set(p1,'Color',colorLines(5,:));
        actClr=sprintf('#%s%s%s',dec2hex(round(colorLines(5,1)*255)),dec2hex(round(colorLines(5,2)*255)),dec2hex(round(colorLines(5,3)*255)));
        results(end+1,:)={colergen(actClr,num2str(meanBOOT(17))),colergen(actClr,num2str(stdBOOT(17))),colergen(actClr,num2str(meanBOOT(18))),colergen(actClr,num2str(stdBOOT(18))),colergen(actClr,num2str(meanBOOT(19)/sqrt(2))),colergen(actClr,num2str(stdBOOT(19)/sqrt(2))),colergen(actClr,num2str(meanBOOT(20))),colergen(actClr,num2str(stdBOOT(20)))};
        results2EXP=[results2EXP;meanBOOT(17),stdBOOT(17),meanBOOT(18),stdBOOT(18),meanBOOT(19)/sqrt(2),stdBOOT(19)/sqrt(2),meanBOOT(20),stdBOOT(20)];        
    end
    
    hTotal=plot(handles.axesFRET,edgesFine+0.5*(edgesFRET(2)-edgesFRET(1)),mean(bootFitAll,2),'--k','LineWidth',2);
    
    if boolTotal
        set(hTotal,'Visible','on');
    else
        set(hTotal,'Visible','off');
    end  
    
    set(handles.uitableResults,'Data',results);
    set(handles.textRMSE,'String',['RMSE: ' num2str(gof.rmse)]);
    set(handles.checkSumFit,'Enable','on');
    set(handles.pushbuttonExport,'Enable','on');
    set(handles.pushbuttonExport,'CData',ButtonLayout2); 
    set(handles.BootSample,'Enable','on');
    set(handles.Probs,'Enable','on');
    set(handles.BootFits,'Enable','on');
    ylim=get(handles.axesFRET,'YLim');
    set(handles.axesFRET,'YLim',[0 ylim(2)]);
    

function checkSumFit_Callback(hObject, eventdata, handles)

    global hTotal;
    
    boolTotal=boolean(get(handles.checkSumFit,'Value'));
    
    if boolTotal
        set(hTotal,'Visible','on');
    else
        set(hTotal,'Visible','off');
    end

function pushbuttonExport_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global weightFRET;
    global sumFRET;
    global edgesFine;
    global cfunX;
    global cfunXMin;
    global cfunXMax;
    global results2EXP;     
    
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    
    if boolWeights
        
        showFRET=weightFRET./(1-boolNorm+boolNorm*sum(weightFRET));
    else
        
        showFRET=sumFRET./(1-boolNorm+boolNorm*sum(sumFRET));
    end
    
    params=get(handles.uitableResults,'Data');

    opt_func=[];
    max_func=[];
    min_func=[];
    
    for iter=1:size(params,1)
        
        opt_func=[opt_func cfunX{iter}];
        min_func=[min_func cfunXMin{iter}];
        max_func=[max_func cfunXMax{iter}];
    end
    
    folder=uigetdir('');
    
    prompt = {'Enter export name:'};
    dlg_title = 'Input';
    num_lines = 1;
    defaultans = {''};
    answer = inputdlg(prompt,dlg_title,num_lines,defaultans);
    
    data.edges=edgesFRET;
    data.histogram=showFRET;
    data.fits.description_results={'Amp','errAmp','Pos','errPos','Sigma','errSigma','Probability','errProbability'};
    data.fits.results=results2EXP;
    data.fits.functions.edges=edgesFine';
    data.fits.functions.opt_func=opt_func;
    data.fits.functions.min_func=min_func;
    data.fits.functions.max_func=max_func;
    
    save([folder '\' answer{1} '.mat'],'data');

function Export_Graph_Callback(hObject, eventdata, handles)

    set(0,'showhiddenhandles','on');
    h=findobj(gcf,'type','axes');
    f1=figure;
    copyobj(h,f1);
    set(gca,'Units','normalized','Position',[0.125 0.1 0.775 0.825],'Box','on');
    set(gcf,'Color',[1 1 1]);

function pushbuttonInput_Callback(hObject, eventdata, handles)

    numGauss=get(handles.popupGauss,'Value')-1;
    [xG,yG]=ginput(numGauss);
    xG=[0;xG;1];
    
    uiData=get(handles.uitableGauss,'Data');
    
    for iter=1:numGauss
        
        uiData(iter,1)=yG(iter);
        uiData(iter,2)=0;
        uiData(iter,3)=sum(yG);
        uiData(iter,4)=xG(iter+1);
        uiData(iter,5)=xG(iter);
        uiData(iter,6)=xG(iter+2);
        uiData(iter,7)=0.1;
        uiData(iter,8)=0;
        uiData(iter,9)=1;
    end
    
    set(handles.uitableGauss,'Data',uiData);

function checkNorm_Callback(hObject, eventdata, handles)

    checkWeight_Callback(hObject, eventdata, handles);
    
    set(handles.BootSample,'Enable','off');
    set(handles.Probs,'Enable','off');
    set(handles.BootFits,'Enable','off');
    set(handles.pushbuttonExport,'CData',[]);    
    set(handles.pushbuttonExport,'Enable','off');
    set(handles.uitableResults,'Data',[]);
    axis 'auto y';
    

function editSamples_Callback(hObject, eventdata, handles)

function editSamples_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function editIterations_Callback(hObject, eventdata, handles)

function editIterations_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function Untitled_1_Callback(hObject, eventdata, handles)

function BootSample_Callback(hObject, eventdata, handles)

    global bootFRET;
    global edgesFRET;
        
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    
    numIter=size(bootFRET,2);
    
    figure();
    hold on;
    
    for iter=1:numIter
        
        plot(edgesFRET,bootFRET(:,iter));
    end
    
    xlim([-0.1 1.1]);
    xlabel('FRET efficiency');
    set(gca,'Box','on');
    
    if boolWeights
        
        if boolNorm
            ylabel('\itp\rm(\Delta\itE\rm)');
        else
            ylabel('Number of molecules');
        end
    else
        if boolNorm
            ylabel('\itp\rm(\Delta\itE\rm)');
        else
            ylabel('Number of frames');
        end
    end    

function BootFits_Callback(hObject, eventdata, handles)

    global edgesFRET;
    global edgesFine;
    global bootFRET;
    global bootFitAll;
        
    boolWeights=boolean(get(handles.checkWeight,'Value'));
    boolNorm=boolean(get(handles.checkNorm,'Value'));
    
    numIter=size(bootFitAll,2);
    
    lmrArray=lrSD(bootFitAll(:,:));
    cfunXMin=lmrArray(:,2)-lmrArray(:,1);
    cfunX=lmrArray(:,2);
    cfunXMax=lmrArray(:,2)+lmrArray(:,3);
    
    lmrArray=lrSD(bootFRET(:,:));
    
    bootXMin=lmrArray(:,2)-lmrArray(:,1);
    bootX=lmrArray(:,2);
    bootXMax=lmrArray(:,2)+lmrArray(:,3);
    
    figure();
    hold on;
    xxArea=[edgesFine(:);flipud(edgesFine(:))];
    yyArea=[cfunXMax(:);flipud(cfunXMin(:))];   
    ff1=fill(xxArea,yyArea,'r');
    set(ff1,'FaceAlpha',0.3,'FaceColor',[1 0 0],'EdgeColor',[1 0 0]);
    plot(edgesFine(:),cfunX,'-r','LineWidth',2);
    errorbar(edgesFRET,bootX,bootX-bootXMin,bootXMax-bootX,'s-k');
    xlim([-0.1 1.1]);
    xlabel('FRET efficiency');
    set(gca,'Box','on');
    legend('conf. interval of fit curves','average fit curve','average bootstrapping sample');
    
    if boolWeights
        
        if boolNorm
            ylabel('\itp\rm(\Delta\itE\rm)');
        else
            ylabel('Number of molecules');
        end
    else
        if boolNorm
            ylabel('\itp\rm(\Delta\itE\rm)');
        else
            ylabel('Number of frames');
        end
    end
    
function lmrArray=lrSD(funArray)

    % funArray = [fun1(:) fun2(:) ... funN(:)]

    numPoints=size(funArray,1);
    
    meanValues=mean(funArray,2);
    
    leftSD=zeros(numPoints,1);
    rightSD=zeros(numPoints,1);

    for iter=1:numPoints
        
        indLeft=squeeze(funArray(iter,:))<meanValues(iter);
        indRight=squeeze(funArray(iter,:))<meanValues(iter);
        
        if ~isempty(indLeft)
            
            leftSD(iter)=sqrt(sum((meanValues(iter)-funArray(iter,indLeft)).^2)/(sum(indLeft)-1));
        else
            leftSD(iter)=0;
        end
            
        if ~isempty(indRight)
        
            rightSD(iter)=sqrt(sum((meanValues(iter)-funArray(iter,indRight)).^2)/(sum(indRight)-1));
        else
            rightSD(iter)=0;
        end
    end
        
%     figure();
%     hold on;
%     
%     for iter=1:size(funArray,2)
%         
%        plot(funArray(:,iter),'-k'); 
%     end
%     
%     plot(meanValues,'-r');
%     hold on;
%     plot(meanValues-leftSD,'--r');
%     plot(meanValues+rightSD,'--r');
    
    lmrArray=[leftSD(:) meanValues(:) rightSD(:)];

function Probs_Callback(hObject, eventdata, handles)

    global results2EXP;
    
    % load colormap for fit curves
    colorLines=colormap('lines');
    
    numGauss=length(results2EXP(:,7));
    
    figure();
    hold on;
    
    for iter=1:numGauss
        
        b1=bar(iter,results2EXP(iter,7));
        set(b1,'FaceColor',colorLines(iter,:));
        errorbar(iter,results2EXP(iter,7),results2EXP(iter,8),'.k');    
    end
    
    xlim([0 numGauss+1]);
    xlabel('State');
    ylabel('Probability');


function figure1_CloseRequestFcn(hObject, eventdata, handles)
    
delete(hObject);
