%% Plot Looking at Category Over Time: Face ROIS
channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
vardir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/';
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';
load([inddir 'SSCM_ROI_Include.mat']);
load([inddir 'SSCM_indOT.mat']);
load([inddir 'SSCM_ROI_All_area.mat']);
load([inddir 'SSCM_Set1.2.3.4.5.6.7_indTrl_ROI_Cat.mat']);
% load([vardir 'SSCM_fruitInd_S1-7.mat']);
clear Sets
load([inddir 'SSCM_ROI_Include.mat']);
load([inddir 'SSCM_indOT.mat']);
%% Plot Probability of Viewing Monkeys & Objects on OT & SL
fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=true;
JN=true;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end

% Load only necessary set of fixHit 
loadSet=1;
if ~loadSet
    load([inddir 'SSCM_fixHit_Whole.mat']);
end

% What Time Period??????
begTimStr=0;
if begTimStr==0
    begTim=1;
    begTimStr=(1/fsamp);
else
    begTim=begTimStr*fsamp;
end
endTimStr=6;
endTim=endTimStr*fsamp;

for pres=1:2;
    for cat=1:size(ROI_Cat{1,1},1)
        D=cell(size(ROI_Cat{1,1}{cat,1},2),2);
        for level=1:size(ROI_Cat{1,1}{cat,1},2);
            ot=1;
            sl=1;
            for setloop=2:7
                if loadSet
                    saveDir=[fixDir num2str(setloop) '/'];
                    load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                end
                for fileloop=1:size(Sets{setloop,1},2)
                    if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                        
                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                        
                        cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                        xMat=cellMat(:,begTim:endTim);
                        clear cellMat
                        
                        xArea=ROI_All_area{setloop,1}(trlInd,1);
                        sceneList=ROI_Scene{setloop,1}(trlInd,1);
                        sceneListU=unique(sceneList);
                        sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                        for k=1:size(sceneListU,1);
                            ind=sceneList==sceneListU(k,1);
                            % Normalize looking time by ROI area if category is one
                            % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
                            
                            if normArea
                                if cat>2 && cat~=6
                                    sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                else
                                    sceneSum(k,:)=mean(xMat(ind,:),1);
                                end
                            else
                                sceneSum(k,:)=mean(xMat(ind,:),1);
                            end
                            
                        end
                        clear xMat
                        
                        if indOT{setloop,1}(1,fileloop);
                            D{level,2}(ot:(ot+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                            ot=ot+(size(sceneSum,1)-1);
                        else
                            D{level,1}(sl:(sl+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                            sl=sl+(size(sceneSum,1)-1);
                        end
                        clear sceneSum
                    end
                end
            end
            if loadSet
                clear fixHit_Whole
            end
        end
    end
end

%% Plot Probability of Viewing Monkeys & Objects on OT & SL, Separate Plot
% for Each Level
fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];

% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;
% Normalize Viewing Time to ROI Area in Unbalanced Conditions?
if normArea
    normStr='Norm';
    yLabel='Normed Probability of Viewing ROI';
elseif normNum
    normStr='Mean';
    yLabel='Probability of Viewing ROI (%)';
end

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

drugLabel={'Saline','Oxytocin'};


% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=0;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end

% Load only necessary set of fixHit 
loadSet=1;
if ~loadSet
    load([inddir 'SSCM_fixHit_Whole.mat']);
end

% Plot SL & OT for each Level Separately
plotDrug=1;

% Plot All Levels Together for OT & SL Separately
plotLevel=1;

% Plot Presentation 1 & 2 Together for Each Level Separately
plotPres=1;

presEnd=[10,6];

roiLabl='Whole';

% Save as .eps without title
saveEPS=0;
savePNG=1;

for pres=1:2;
    
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    
    for cat=1:size(ROI_Cat{1,1},1)
        if normArea && cat==1
            yLabel=['Normed ',yLabel];
        elseif normArea && cat~=1
            yLabel='Probability of Viewing ROI (%)';
        end
        D=cell(size(ROI_Cat{1,1}{cat,1},2),2);
        for level=1:size(ROI_Cat{1,1}{cat,1},2);
            ot=1;
            sl=1;
            for setloop=2:7
                if loadSet
                    saveDir=[fixDir num2str(setloop) '/'];
                    load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                end
                for fileloop=1:size(Sets{setloop,1},2)
                    if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                        
                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                        
                        cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                        xMat=cellMat(:,begTim:endTim);
                        clear cellMat
                        
                        xArea=ROI_All_area{setloop,1}(trlInd,1);
                        sceneList=ROI_Scene{setloop,1}(trlInd,1);
                        sceneListU=unique(sceneList);
                        sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                        for k=1:size(sceneListU,1);
                            ind=sceneList==sceneListU(k,1);
                            % Normalize looking time by ROI area if category is one
                            % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
                            
                            if normArea
                                if cat>2 && cat~=6
                                    sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                else
                                    sceneSum(k,:)=mean(xMat(ind,:),1);
                                end
                            elseif normNum
                                if cat==1
                                    sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                else
                                    sceneSum(k,:)=mean(xMat(ind,:),1);
                                end
                            end
                            
                        end
                        clear xMat
                        
                        if indOT{setloop,1}(1,fileloop);
                            D{level,2}(ot:(ot+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                            ot=ot+size(sceneSum,1);
                        else
                            D{level,1}(sl:(sl+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                            sl=sl+size(sceneSum,1);
                        end
                        clear sceneSum
                    end
                end
            end
        end 
        
        if loadSet
            clear fixHit_Whole
        end
        clear channelcolormap
        channelcolormap{1} = [.2 .8 .2;0 .5 0];
        channelcolormap{2} = [1 0 0;.545 0 0];
        channelcolormap{3} = [.11 .56 1;0 0 .545];
%         channelcolormap = [1 .55 0;1 .27 0];

        % Plot SL & OT for each Level Separately
        if plotDrug
            yMax=0;
            for drug=1:2
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    maxCur=(max(mean(D{level,drug},1))*100)*1.1;
                    if maxCur>yMax
                        yMax=round(maxCur);
                    end
                end
            end
            for level=1:size(ROI_Cat{1,1}{cat,1},2);
                figure;
%                 for drug=1:2
%                 dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},channelcolormap(drug,:),1,smval,0,0,1,fsamp,0)
%                 hold on;
%                 n{drug}=[drugLabel{drug},', N=',num2str(size(D{level,drug},1))];
%                 end
                for drug=1:2
                    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},channelcolormap{level}(drug,:),1,smval,0,0,1,fsamp,0)
                    hold on;
                    n{drug}=[drugLabel{drug},', N=',num2str(size(D{level,drug},1))];
                end
                titleLabel=['Whole Item, ',ROI_Cat_Labl{cat,1},', ',ROI_Cat_Labl{cat,2}{1,level},', ',presLabel{pres},', ',normStr];
                title(titleLabel);
                legend(n,'Location','NorthEast');
%                 yLim([0 yMax]);
                yLim([0 25]);
                ylabel(yLabel);
                xlabel('Time From Stimulus Onset (s)');
                
                if MP && ~JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif JN && ~MP
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif MP && JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                else
                end
                
                if ~exist(newDir,'dir');
                    mkdir(newDir);
                end
                
                saveName=[newDir,'PSTH_Set2-7_OT-SL,',mStr,'_',titleLabel];
                
                if saveEPS
                    if savePNG
                        export_fig([saveName,'.eps'],'-eps')
                        title(titleLabel)
                        export_fig([saveName,'.png'],'-png')
                    end
                elseif savePNG
                    title(titleLabel)
                    export_fig([saveName,'.png'],'-png')
                else
                end
                close(gcf)
            end
        end
        
        % Plot All Levels Together for OT & SL Separately
        if plotLevel
            yMax=0;
            for drug=1:2
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    maxCur=(max(mean(D{level,drug},1))*100)*1.05;
                    if maxCur>yMax
                        yMax=round(maxCur);
                    end
                end
            end
            for drug=1:2;
                figure;
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},channelcolormap{level}(drug,:),1,smval,0,0,1,fsamp,0)                    
                    hold on;
                    n{level}=[ROI_Cat_Labl{cat,2}{1,level},', N=',num2str(size(D{level,drug},1))];
                end
                titleLabel=['Whole Item, ',ROI_Cat_Labl{cat,1},', ',drugLabel{drug},', ',presLabel{pres},', ',normStr];
                title(titleLabel);
                
                legend(n,'Location','NorthEast');
                yLim([0 yMax]);
                ylabel(yLabel);
                xlabel('Time From Stimulus Onset (s)');
                
                if MP && ~JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif JN && ~MP
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif MP && JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                else
                end
                
                if ~exist(newDir,'dir');
                    mkdir(newDir);
                end
                
                saveName=[newDir,'PSTH_Set2-7_OT-SL,',mStr,'_',titleLabel];
                
                if saveEPS
                    if savePNG
                        export_fig([saveName,'.eps'],'-eps')
                        title(titleLabel)
                        export_fig([saveName,'.png'],'-png')
                    end
                elseif savePNG
                    title(titleLabel)
                    export_fig([saveName,'.png'],'-png')
                else
                end
                close(gcf)
            end
        end        
    end
end

%% Plot Probability of Viewing Monkeys (Age) & Objects on OT & SL, Separate Plot
% for Each Level, FOR Age+Object RASTER
%% Raster Plot for Monkeys & Objects, pool subjects
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% # of scenes in Set
nScenes=90;
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=1;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);
% Area of ROI
% channelcolormap = [.25 .4 .88;.6 0 0;0 .4 0;1 .6 0];
% channelcolormap = [.11 .56 1;.9 0 0;.2 .8 .2;1 .27 0];
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];


presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=3;

% List of categories to plot
catList=[4,4,4,2];
levelList=[3,2,1,2];

sortStr=ROI_Cat_Labl{catList(sortBy),2}{levelList(sortBy)};

% Save as .eps without title
saveEPS=1;
savePNG=1;



for pres=1:2;
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
        
    
    for drug=1:2        
        D=cell(length(levelList),2);
        c=1;
        shift=0;
        for setloop=2:7
            saveDir=[fixDir num2str(setloop) '/'];
            load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
            for fileloop=1:size(Sets{setloop,1},2)
                if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                    if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                        
                            for levelLoop=1:length(levelList)
                                level=levelList(levelLoop);
                                cat=catList(levelLoop);                                
                                
                                trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                
                                cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                xMat=cellMat(:,begTim:endTim);
                                clear cellMat
                                
                                xArea=ROI_All_area{setloop,1}(trlInd,1);
                                sceneList=ROI_Scene{setloop,1}(trlInd,1);
                                sceneListU=unique(sceneList);
                                sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                                for k=1:size(sceneListU,1);
                                    ind=sceneList==sceneListU(k,1);
                                    % Normalize looking time by ROI area if category is one
                                    % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
                                    
                                    if normArea
                                        if cat>2 && cat~=6
                                            sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                        else
                                            sceneSum(k,:)=mean(xMat(ind,:),1);
                                        end
                                    elseif normNum
                                        if cat==1
                                            sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                        else
                                            sceneSum(k,:)=mean(xMat(ind,:),1);
                                        end
                                    end
                                    
                                end
                                clear xMat
                                
                                
                                if indOT{setloop,1}(1,fileloop);
                                    D{levelLoop,drug}(c,begTim:endTim)=sceneSum;
                                else
                                    D{levelLoop,drug}(c,begTim:endTim)=zeros(1,endTim);
                                end
                                clear sceneSum
                            
                        end
                    end
                end
            end
            clear fixHit_Whole
        end
        % Sort trials by time spent in level 1 ROI
        [~,sortInd]=sort(sum(D{sortBy,drug},2),1,'descend');
        figure;
        for levelLoop=1:length(levelList)
            raster2(D{levelLoop,drug}(sortInd,:),channelcolormap(levelLoop,:),shift);
            hold on;
        end
        titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug},'Sorted By ',sortStr];
        ylim([0 size(D{1,drug},1)]);
        set(gca,'YDir','reverse');
        clear D

        %         B = axes;
        %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
        %         ylabel(B,'Right Label')
        newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
        saveName=[newdir,'ROI_Set2-7_',titleLabel];
        
        if ~exist(newdir,'dir');
            mkdir(newdir);
        end
        
        if saveEPS
            if savePNG
                export_fig([newdir,'.eps'],'-eps')
                saveas(gcf,[saveName,'.eps'],'eps')
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            end
        elseif savePNG
            title(titleLabel)
            saveas(gcf,[saveName,'.png'],'png')
        else
        end
        close(gcf)
    end
end
%%

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];

% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;
% Normalize Viewing Time to ROI Area in Unbalanced Conditions?
if normArea
    normStr='Norm';
    yLabel='Normed Probability of Viewing ROI';
elseif normNum
    normStr='Mean';
    yLabel='Probability of Viewing ROI (%)';
end

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

drugLabel={'Saline','Oxytocin'};


% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=true;
JN=true;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end

% Load only necessary set of fixHit 
loadSet=1;
if ~loadSet
    load([inddir 'SSCM_fixHit_Whole.mat']);
end

% Plot SL & OT for each Level Separately
plotDrug=1;

% Plot All Levels Together for OT & SL Separately
plotLevel=1;

% Plot Presentation 1 & 2 Together for Each Level Separately
plotPres=1;

presEnd=[10,6];

roiLabl='Whole';


% List of categories to plot
catList=[4,4,4,2];
levelList=[3,2,1,2];

% Save as .eps without title
saveEPS=1;
savePNG=1;

for pres=2;
    
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    
 
        D=cell(length(levelLoop),2);
            
            for setloop=2:7
                if loadSet
                    saveDir=[fixDir num2str(setloop) '/'];
                    load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                end
                for fileloop=1:size(Sets{setloop,1},2)
                    if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                        for levelLoop=1:length(levelList)
                            level=levelList(levelLoop);
                            cat=catList(levelLoop);
                            
                            trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                            
                            cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                            xMat=cellMat(:,begTim:endTim);
                            clear cellMat
                            
                            xArea=ROI_All_area{setloop,1}(trlInd,1);
                            sceneList=ROI_Scene{setloop,1}(trlInd,1);
                            sceneListU=unique(sceneList);
                            sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                            for k=1:size(sceneListU,1);
                                ind=sceneList==sceneListU(k,1);
                                % Normalize looking time by ROI area if category is one
                                % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
                                
                                if normArea
                                    if cat>2 && cat~=6
                                        sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                    else
                                        sceneSum(k,:)=mean(xMat(ind,:),1);
                                    end
                                elseif normNum
                                    if cat==1
                                        sceneSum(k,:)=((mean(xMat(ind,:),1))./(sum(xArea(ind,:))/totalArea));
                                    else
                                        sceneSum(k,:)=mean(xMat(ind,:),1);
                                    end
                                end
                                
                            end
                            clear xMat
                            
                            if indOT{setloop,1}(1,fileloop);
                                D{levelLoop,2}(ot:(ot+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                                ot=ot+size(sceneSum,1);
                            else
                                D{levelLoop,1}(sl:(sl+(size(sceneSum,1)-1)),begTim:endTim)=sceneSum;
                                sl=sl+size(sceneSum,1);
                            end
                            clear sceneSum
                        end
                    end
                end
            end
        end 
        
        if loadSet
            clear fixHit_Whole
        end
        
        
        % Plot SL & OT for each Level Separately
        if plotDrug
            yMax=0;
            for drug=1:2
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    maxCur=(max(mean(D{level,drug},1))*100)*1.1;
                    if maxCur>yMax
                        yMax=round(maxCur);
                    end
                end
            end
            for level=1:size(ROI_Cat{1,1}{cat,1},2);
                figure;
                for drug=1:2
                dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},channelcolormap(drug,:),1,smval,0,0,1,fsamp,0)
                hold on;
                n{drug}=[drugLabel{drug},', N=',num2str(size(D{level,drug},1))];
                end
                titleLabel=['Whole Item, ',ROI_Cat_Labl{cat,1},', ',ROI_Cat_Labl{cat,2}{1,level},', ',presLabel{pres},', ',normStr];
                title(titleLabel);
                legend(n,'Location','NorthEast');
                yLim([0 yMax]);
                ylabel(yStr);
                xlabel('Time From Stimulus Onset (s)');
                
                if MP && ~JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif JN && ~MP
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif MP && JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                else
                end
                
                if ~exist(newDir,'dir');
                    mkdir(newDir);
                end
                
                saveName=[newDir,'PSTH_Set2-7_OT-SL,',mStr,'_',titleLabel];
                
                if saveEPS
                    if savePNG
                        export_fig([saveName,'.eps'],'-eps')
                        title(titleLabel)
                        export_fig([saveName,'.png'],'-png')
                    end
                elseif savePNG
                    title(titleLabel)
                    export_fig([saveName,'.png'],'-png')
                else
                end
                close(gcf)
            end
        end
        
        % Plot All Levels Together for OT & SL Separately
        if plotLevel
            yMax=0;
            for drug=1:2
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    maxCur=(max(mean(D{level,drug},1))*100)*1.05;
                    if maxCur>yMax
                        yMax=round(maxCur);
                    end
                end
            end
            for drug=1:2;
                figure;
                for level=1:size(ROI_Cat{1,1}{cat,1},2);
                    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},channelcolormap(level,:),1,smval,0,0,1,fsamp,0)
                    hold on;
                    n{level}=[ROI_Cat_Labl{cat,2}{1,level},', N=',num2str(size(D{level,drug},1))];
                end
                titleLabel=['Whole Item, ',ROI_Cat_Labl{cat,1},', ',drugLabel{drug},', ',presLabel{pres},', ',normStr];
                title(titleLabel);
                
                legend(n,'Location','NorthEast');
                yLim([0 yMax]);
                ylabel(yStr);
                xlabel('Time From Stimulus Onset (s)');
                
                if MP && ~JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif JN && ~MP
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                elseif MP && JN
                    newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
                else
                end
                
                if ~exist(newDir,'dir');
                    mkdir(newDir);
                end
                
                saveName=[newDir,'PSTH_Set2-7_OT-SL,',mStr,'_',titleLabel];
                
                if saveEPS
                    if savePNG
                        export_fig([saveName,'.eps'],'-eps')
                        title(titleLabel)
                        export_fig([saveName,'.png'],'-png')
                    end
                elseif savePNG
                    title(titleLabel)
                    export_fig([saveName,'.png'],'-png')
                else
                end
                close(gcf)
            end
        end        
    end
end


%% Raster Plot for Monkeys & Objects
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=0;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);

presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=[1,1,nan,nan,nan,1];

% List of categories to plot
catList=[1,2,6];

% Save as .eps without title
saveEPS=1;
savePNG=1;


for catLoop=1:size(catList,2);
    cat=catList(catLoop);
    for pres=1:2;
        % What Time Period??????
        begTimStr=0;
        if begTimStr==0
            begTim=1;
            begTimStr=(1/fsamp);
        else
            begTim=begTimStr*fsamp;
        end
        endTimStr=presEnd(pres);
        endTim=endTimStr*fsamp;
        
        for drug=1:2
            figure;
            shift=0;
            lim=20;
            for m=1:size(mList,2);
                mStr=mList{m};
                D=cell(size(ROI_Cat{1,1}{cat,1},2),2);
                for level=1:size(ROI_Cat{1,1}{cat,1},2)
                    c=1;
                    for setloop=2:7
                        saveDir=[fixDir num2str(setloop) '/'];
                        load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                        for fileloop=1:size(Sets{setloop,1},2)
                            if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                                if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                                    % When analyzing ROI Area, only look at monkeys
                                    if cat==1
                                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    else
                                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    end
                                    cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                    xMat=cellMat(:,begTim:endTim);
                                    clear cellMat
                                    
                                    xArea=ROI_All_area{setloop,1}(trlInd,1);
                                    sceneList=ROI_Scene{setloop,1}(trlInd,1);
                                    sceneListU=unique(sceneList);
                                    sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                                    
                                    for k=1:size(sceneListU,1);
                                        ind=sceneList==sceneListU(k,1);
                                        sceneSum(k,:)=sum(xMat(ind,:),1);
                                    end
                                    clear xArea
                                    clear xMat
                                    D{level,drug}(c:(c+(size(sceneSum,1)-1)),:)=sceneSum;
                                    
                                    c=c+size(sceneSum,1);
                                    clear sceneSum
                                end
                            end
                        end
                        clear fixHit_Whole
                    end
                    % Sort trials by time spent in level 1 ROI
                    [~,sortInd]=sort(sum(D{sortBy(cat),drug},2),1,'descend');
                    raster2(D{level,drug}(sortInd,:),channelcolormap(level,:),shift);
                    hold on;
                end
                clear D
                shift=c;
                if m<size(mList,2)
                    shift=shift+5;
                    line([0 endTim],[shift shift],'color','k','linewidth',2,'linestyle','--')
                    shift=shift+7;
                end
                lim=lim+c;
            end
            
            titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug}];
            ylim([0 lim]);
            set(gca,'YDir','reverse');
            %         B = axes;
            %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
            %         ylabel(B,'Right Label')
            newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
            saveName=[newdir,'ROI_Set2-7_',titleLabel];
            
            if ~exist(newdir,'dir');
                mkdir(newdir);
            end
            
            if saveEPS
                if savePNG
                    saveas(gcf,[saveName,'.eps'],'eps')
                    title(titleLabel)
                    saveas(gcf,[saveName,'.png'],'png')
                end                
            elseif savePNG
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            else
            end
            close(gcf)
        end
    end
end

%% Raster Plot for Monkeys & Objects, 
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=0;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);

presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=[1,1,nan,nan,nan,1];

% List of categories to plot
catList=[2,4];
% List of levels to plot
levelList{1}=2;
levelList{2}=[1,2,3];
numLevels=0;
for k=1:length(levelList)
    numLevels=numLevels+numel(levelList{k});
end
% Save as .eps without title
saveEPS=1;
savePNG=1;



for pres=1:2;
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    
    for drug=1:2
        figure;
        shift=0;
        lim=20;
        for m=1:size(mList,2);
            mStr=mList{m};
            D=cell(numLevels,2);
            for catLoop=1:size(catList,2);
                cat=catList(catLoop);
                for levelLoop=1:length(levelList{catLoop})
                    level=levelList{catLoop}(levelLoop);
                    c=1;
                    for setloop=2:7
                        saveDir=[fixDir num2str(setloop) '/'];
                        load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                        for fileloop=1:size(Sets{setloop,1},2)
                            if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                                if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                                    % When analyzing ROI Area, only look at monkeys
                                    if cat==1
                                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    else
                                        trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    end
                                    cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                    xMat=cellMat(:,begTim:endTim);
                                    clear cellMat
                                    
                                    xArea=ROI_All_area{setloop,1}(trlInd,1);
                                    sceneList=ROI_Scene{setloop,1}(trlInd,1);
                                    sceneListU=unique(sceneList);
                                    sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                                    
                                    for k=1:size(sceneListU,1);
                                        ind=sceneList==sceneListU(k,1);
                                        sceneSum(k,:)=sum(xMat(ind,:),1);
                                    end
                                    clear xArea
                                    clear xMat
                                    D{level,drug}(c:(c+(size(sceneSum,1)-1)),:)=sceneSum;
                                    
                                    c=c+size(sceneSum,1);
                                    clear sceneSum
                                end
                            end
                        end
                        clear fixHit_Whole
                    end
                    % Sort trials by time spent in level 1 ROI
                    [~,sortInd]=sort(sum(D{sortBy(cat),drug},2),1,'descend');
                    raster2(D{level,drug}(sortInd,:),channelcolormap(level,:),shift);
                    hold on;
                end
                clear D
                shift=c;
                if m<size(mList,2)
                    shift=shift+5;
                    line([0 endTim],[shift shift],'color','k','linewidth',2,'linestyle','--')
                    shift=shift+7;
                end
                lim=lim+c;
            end
            
            titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug}];
            ylim([0 lim]);
            set(gca,'YDir','reverse');
            %         B = axes;
            %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
            %         ylabel(B,'Right Label')
            newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
            saveName=[newdir,'ROI_Set2-7_',titleLabel];
            
            if ~exist(newdir,'dir');
                mkdir(newdir);
            end
            
            if saveEPS
                if savePNG
                    saveas(gcf,[saveName,'.eps'],'eps')
                    title(titleLabel)
                    saveas(gcf,[saveName,'.png'],'png')
                end
            elseif savePNG
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            else
            end
            close(gcf)
        end
    end
end
%% Raster Plot for Monkeys & Objects, separate by subject
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% # of scenes in Set
nScenes=90;
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=1;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);
% Area of ROI
% channelcolormap = [.25 .4 .88;.6 0 0;0 .4 0;1 .6 0];
% channelcolormap = [.11 .56 1;.9 0 0;.2 .8 .2;1 .27 0];
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];


presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=3;

% List of categories to plot
catList=[4,4,4,2];
levelList=[3,2,1,2];

% Save as .eps without title
saveEPS=1;
savePNG=1;



for pres=1:2;
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    
    for drug=1:2
        figure;
        shift=0;
        lim=20;
        for m=1:size(mList,2);
            mStr=mList{m};
            D=cell(length(levelList),2);
            c=1;
            for setloop=2:7
                saveDir=[fixDir num2str(setloop) '/'];
                load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                for fileloop=1:size(Sets{setloop,1},2)
                    if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                        if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                            for s=1:nScenes
                                sceneInd=ROI_Scene{setloop,1}==s;
                                for levelLoop=1:length(levelList)
                                    level=levelList(levelLoop);
                                    cat=catList(levelLoop);
                                    % When analyzing ROI Area, only look at monkeys
                                    if cat==1
                                        trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    else
                                        trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                    end
                                    if sum(trlInd)>0
                                        cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                        xMat=cellMat(:,begTim:endTim);
                                        clear cellMat
                                        sceneSum=sum(xMat,1);
                                        clear xMat
                                        D{levelLoop,drug}(c,begTim:endTim)=sceneSum;
                                    else
                                        D{levelLoop,drug}(c,begTim:endTim)=zeros(1,endTim);
                                    end
                                    clear sceneSum
                                end
                                c=c+1;
                            end
                        end
                    end
                end
                clear fixHit_Whole
            end
            % Sort trials by time spent in level 1 ROI
            [~,sortInd]=sort(sum(D{sortBy,drug},2),1,'descend');
            for levelLoop=1:length(levelList)
                raster2(D{levelLoop,drug}(sortInd,:),channelcolormap(levelLoop,:),shift);
                hold on;
            end
            clear D
            shift=c;
            if m<size(mList,2)
                shift=shift+5;
                line([0 endTim],[shift shift],'color','k','linewidth',2,'linestyle','--')
                shift=shift+7;
            end
            lim=lim+c;
        end
        
        titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug}];
        ylim([0 lim]);
        set(gca,'YDir','reverse');
        %         B = axes;
        %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
        %         ylabel(B,'Right Label')
        newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
        saveName=[newdir,'ROI_Set2-7_',titleLabel];
        
        if ~exist(newdir,'dir');
            mkdir(newdir);
        end
        
        if saveEPS
            if savePNG
                saveas(gcf,[saveName,'.eps'],'eps')
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            end
        elseif savePNG
            title(titleLabel)
            saveas(gcf,[saveName,'.png'],'png')
        else
        end
        close(gcf)
    end
end

%% Raster Plot for Monkeys & Objects, pool subjects
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% # of scenes in Set
nScenes=90;
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=1;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);
% Area of ROI
% channelcolormap = [.25 .4 .88;.6 0 0;0 .4 0;1 .6 0];
% channelcolormap = [.11 .56 1;.9 0 0;.2 .8 .2;1 .27 0];
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];


presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=3;

% List of categories to plot
catList=[4,4,4,2];
levelList=[3,2,1,2];

sortStr=ROI_Cat_Labl{catList(sortBy),2}{levelList(sortBy)};


% Save as .eps without title
saveEPS=1;
savePNG=1;



for pres=1:2;
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    
    for drug=1:2        
        shift=0;
        D=cell(length(levelList),2);
        trlIndAll=cell(1,length(levelList));
        c=1;
        for setloop=2:7
            saveDir=[fixDir num2str(setloop) '/'];
            load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
            for fileloop=1:size(Sets{setloop,1},2)
                if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                    if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                        for s=1:nScenes
                            sceneInd=ROI_Scene{setloop,1}==s;
                            for levelLoop=1:length(levelList)
                                level=levelList(levelLoop);
                                cat=catList(levelLoop);
                                % When analyzing ROI Area, only look at monkeys
                                if cat==1
                                    trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                else
                                    trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                end
                                trlIndAll{levelLoop}(c,1)=sum(trlInd);
                                if sum(trlInd)>0
                                    cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                    xMat=cellMat(:,begTim:endTim);
                                    clear cellMat
                                    sceneSum=sum(xMat,1);
                                    clear xMat
                                    D{levelLoop,drug}(c,begTim:endTim)=sceneSum;
                                else
                                    D{levelLoop,drug}(c,begTim:endTim)=zeros(1,endTim);
                                end
                                clear sceneSum
                            end
                            c=c+1;
                        end
                    end
                end
            end
            clear fixHit_Whole
        end
        sortSum=sum(D{sortBy,drug},2);
        sortTrls=sortSum>0;
        % Sort trials by time spent in level 1 ROI
        [~,sortInd]=sort(sum(D{sortBy,drug},2),1,'descend');
        [~,sortInd2]=sort(sum(D{1,drug},2),1,'descend');
        figure;
        for levelLoop=1:length(levelList)
            raster2(D{levelLoop,drug}(sortInd,:),channelcolormap(levelLoop,:),shift);
            hold on;
        end
        titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug},'Sorted by ',sortStr];
        ylim([0 size(D{1,drug},1)]);
        set(gca,'YDir','reverse');
        clear D

        %         B = axes;
        %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
        %         ylabel(B,'Right Label')
        newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
        saveName=[newdir,'ROI_Set2-7_',titleLabel];
        
        if ~exist(newdir,'dir');
            mkdir(newdir);
        end
        
        if saveEPS
            if savePNG
                saveas(gcf,[saveName,'.eps'],'eps')
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            end
        elseif savePNG
            title(titleLabel)
            saveas(gcf,[saveName,'.png'],'png')
        else
        end
        close(gcf)
    end
end

%% Raster Plot for Monkeys & Objects, pool subjects, plot levels separately
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% # of scenes in Set
nScenes=90;
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=1;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);
% Area of ROI
% channelcolormap = [.25 .4 .88;.6 0 0;0 .4 0;1 .6 0];
% channelcolormap = [.11 .56 1;.9 0 0;.2 .8 .2;1 .27 0];
clear channelcolormap
channelcolormap = [.11 .56 1;1 0 0;.2 .8 .2;1 .55 0];
%             dodger blue    red,  green, orange

presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=3;

% List of categories to plot
catList=[4,4,4,2];
levelList=[3,2,1,2];
% Adult, Juvenile, Infant, Object
sortStr=ROI_Cat_Labl{catList(sortBy),2}{levelList(sortBy)};


% Save as .eps without title
saveEPS=1;
savePNG=1;



for pres=1:2;
    % What Time Period??????
    begTimStr=0;
    if begTimStr==0
        begTim=1;
        begTimStr=(1/fsamp);
    else
        begTim=begTimStr*fsamp;
    end
    endTimStr=presEnd(pres);
    endTim=endTimStr*fsamp;
    D=cell(length(levelList),2);
    trlIndAll=cell(length(levelList),2);
    for drug=1:2        
        shift=0;
        c=1;
        for setloop=2:7
            saveDir=[fixDir num2str(setloop) '/'];
            load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
            for fileloop=1:size(Sets{setloop,1},2)
                if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                    if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                        for s=1:nScenes
                            sceneInd=ROI_Scene{setloop,1}==s;
                            for levelLoop=1:length(levelList)
                                level=levelList(levelLoop);
                                cat=catList(levelLoop);
                                % When analyzing ROI Area, only look at monkeys
                                if cat==1
                                    trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                else
                                    trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                end
                                trlIndAll{levelLoop,drug}(c,1)=sum(trlInd);
                                if sum(trlInd)>0
                                    cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                    xMat=cellMat(:,begTim:endTim);
                                    clear cellMat
                                    sceneSum=sum(xMat,1);
                                    clear xMat
                                    D{levelLoop,drug}(c,begTim:endTim)=sceneSum;
                                else
                                    D{levelLoop,drug}(c,begTim:endTim)=zeros(1,endTim);
                                end                                
                                clear sceneSum
                            end
                            c=c+1;
                        end
                    end
                end
            end
            clear fixHit_Whole
        end
    end
        
        for levels=1:length(levelList)   
            level=levelList(levels);
            cat=catList(levels);
            figure;
            shift=0;
            for drug=1:2
                hasROI=find(trlIndAll{levels,drug}>0);
                % Sort trials by time spent in level 1 ROI
                [~,sortInd]=sort(sum(D{levels,drug},2),1,'descend');
                sortInd=sortInd(ismember(sortInd,hasROI),1);
                
                for levelLoop=1:length(levelList)
                    raster2(D{levelLoop,drug}(sortInd,:),channelcolormap(levelLoop,:),shift);
                    hold on;
                end
                shift=shift+size(sortInd,1);
                if drug==1
                    shift=shift+round(size(sortInd,1)*.02);
                    line([0 endTim],[shift shift],'color','k','linewidth',2,'linestyle','--')
                    shift=shift+round(size(sortInd,1)*.022);
                end
            end
            
            titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,2}{level},', ',presLabel{pres},', MP & JN, SL & OT'];
            ylim([0 shift]);
            set(gca,'YDir','reverse');
            newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
            saveName=[newdir,'ROI_Set2-7_',titleLabel];
            
            if ~exist(newdir,'dir');
                mkdir(newdir);
            end
            
            if saveEPS
                if savePNG
                    saveas(gcf,[saveName,'.eps'],'eps')
                    title(titleLabel)
                    saveas(gcf,[saveName,'.png'],'png')
                end
            elseif savePNG
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            else
            end
            close(gcf)
        end 
%         
%         % Plot SL & OT for each Level Separately
%         if plotDrug
%             yMax=0;
%             for drug=1:2
%                 for levels=1:length(levelList)
%                     maxCur=(max(mean(D{levels,drug},1))*100)*1.1;
%                     if maxCur>yMax
%                         yMax=round(maxCur);
%                     end
%                 end
%             end
%             for levels=1:length(levelList)   
%                 level=levelList(levels);
%                 cat=catList(levels);
%                 figure;
%                 for drug=1:2
%                     % Sort trials by time spent in level 1 ROI
%                     [~,sortInd]=sort(sum(D{levels,drug},2),1,'descend');
%                     sortIndTrls=sortInd(trlIndAll{levels}>0);
% %                     dofill((begTimStr:(1/fsamp):endTimStr),D{levels,drug}(sortIndTrls,:),channelcolormap(levels,:),1,smval,0,0,1,fsamp,0)
%                     dofill((begTimStr:(1/fsamp):endTimStr),D{levels,drug},channelcolormap(levels,:),1,smval,0,0,1,fsamp,0)
%                     hold on;
%                     n{drug}=[drugLabel{drug},', N=',num2str(sum(trlIndAll{levels}>0))];
%                 end
% %                 titleLabel=['Whole Item, ',ROI_Cat_Labl{cat,1},', ',ROI_Cat_Labl{cat,2}{1,level},', ',presLabel{pres},', ',normStr];
% %                 title(titleLabel);
% %                 legend(n,'Location','NorthEast');
% % %                 yLim([0 yMax]);
% %                 yLim([0 25]);
% %                 ylabel(yLabel);
% %                 xlabel('Time From Stimulus Onset (s)');
%                 
% %                 if MP && ~JN
% %                     newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
% %                 elseif JN && ~MP
% %                     newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
% %                 elseif MP && JN
% %                     newDir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/',roiLabl,'ROIs/',presLabel{pres},'/Fixation/',yLabel,'/'];
% %                 else
% %                 end
% %                 
% %                 if ~exist(newDir,'dir');
% %                     mkdir(newDir);
% %                 end
% %                 
% %                 saveName=[newDir,'PSTH_Set2-7_OT-SL,',mStr,'_',titleLabel];
% %                 
% %                 if saveEPS
% %                     if savePNG
% %                         export_fig([saveName,'.eps'],'-eps')
% %                         title(titleLabel)
% %                         export_fig([saveName,'.png'],'-png')
% %                     end
% %                 elseif savePNG
% %                     title(titleLabel)
% %                     export_fig([saveName,'.png'],'-png')
% %                 else
% %                 end
% %                 close(gcf)
%             end
%         end
        clear D    
end
%% Raster Plot for Monkeys & Objects
clear fruitInd trlType trlType_Labl

fixDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/fixHit_Whole';
% # of scenes in Set
nScenes=90;
% Novel Only
smval=200;
%Sampling Frequency
fsamp=1000;
%Total Area of Image in Pixels
totalArea=800*600;
% Normalize by Cumulative % of Image Occupied by All ROIs
normArea=0;
% Normalize by # of ROIs (mean)
normNum=1;

% 1st Presentation: Pres=1
% 2nd Presentation: Pres=2;
presLabel=[{'Pres1'},{'Pres2'}];

% Only Objects: Fruit??????
fruit=true;

% MP, JN or MP & JN Combined??????
MP=0;
JN=1;
if MP && ~JN
    mStr='MP';
elseif JN && ~MP
    mStr='JN';
elseif MP && JN
    mStr='MPJN';
else
end
mList={'JN','MP'};

% channelcolormap = [0 0 1;0.75 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap = [0 0 1;1 0 0;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap=nicejet(3);
% channelcolormap=winter(3);
% Area of ROI
channelcolormap{1} = [0 0 1;1 0 0;0 1 0;];
% Item Type
channelcolormap{2} = [0 0 1;1 0 0;];
% # of Eyes Visible
channelcolormap{3} = [0 1 0;1 0 0;0 0 1;];
% Age
channelcolormap{4} = [0 1 0;1 0 0;0 0 1;];
% Sex
channelcolormap{5} = [1 0 0;0 0 1;0 1 0;];
% Gaze Direction
channelcolormap{6} = [0 0 1;1 0 0;];

presEnd=[10,6];

drugLabel={'SL','OT'};
drugNum=[0,1];

% What level should the rasters be sorted by?
sortBy=[1,1,3,1,2,1];

% List of categories to plot
catList=[1,2,3,4,5,6];

% Save as .eps without title
saveEPS=1;
savePNG=1;


for catLoop=5:size(catList,2);
    cat=catList(catLoop);
    if cat==5
        levelEnd=2;
    else
        levelEnd=size(ROI_Cat{1,1}{cat,1},2);
    end
    for pres=1:2;
        % What Time Period??????
        begTimStr=0;
        if begTimStr==0
            begTim=1;
            begTimStr=(1/fsamp);
        else
            begTim=begTimStr*fsamp;
        end
        endTimStr=presEnd(pres);
        endTim=endTimStr*fsamp;
        
        for drug=1:2
            figure;
            shift=0;
            lim=20;
            for m=1:size(mList,2);
                mStr=mList{m};
                D=cell(size(ROI_Cat{1,1}{cat,1},2),2);    
                c=1;
                for setloop=2:7
                    saveDir=[fixDir num2str(setloop) '/'];
                    load([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat']);
                    for fileloop=1:size(Sets{setloop,1},2)
                        if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                            if indOT{setloop,1}(1,fileloop)==drugNum(drug);
                                for s=1:nScenes
                                    sceneInd=ROI_Scene{setloop,1}==s;
                                    for level=1:levelEnd
                                        % When analyzing ROI Area, only look at monkeys
                                        if cat==1
                                            trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & ROI_Cat{setloop,1}{2,1}(:,1) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                        else
                                            trlInd=sceneInd & ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,pres) & ~indTrl{setloop,1}(:,8) & ~indTrl{setloop,1}(:,9) & include{setloop,1}(:,fileloop);
                                        end
                                        if sum(trlInd)>0
                                            cellMat=cell2mat(fixHit_Whole{setloop,1}{1,fileloop}(trlInd,pres));
                                            xMat=cellMat(:,begTim:endTim);
                                            clear cellMat
                                            sceneSum=sum(xMat,1);                                            
                                            clear xMat
                                            D{level,drug}(c,begTim:endTim)=sceneSum;                                             
                                        else
                                            D{level,drug}(c,begTim:endTim)=zeros(1,endTim);
                                        end
                                        clear sceneSum
                                    end
                                    c=c+1;
                                end
                            end
                        end
                    end
                    clear fixHit_Whole
                end
                % Sort trials by time spent in level 1 ROI
                [~,sortInd]=sort(sum(D{sortBy(cat),drug},2),1,'descend');
                for level=1:levelEnd               
                    raster2(D{level,drug}(sortInd,:),channelcolormap{cat}(level,:),shift);
                    hold on;
                end
                clear D
                shift=c;
                if m<size(mList,2)
                    shift=shift+5;
                    line([0 endTim],[shift shift],'color','k','linewidth',2,'linestyle','--')
                    shift=shift+7;
                end
                lim=lim+c;
            end
            
            titleLabel=['Fixation Raster: ',ROI_Cat_Labl{cat,1},', ',presLabel{pres},', MP & JN, ',drugLabel{drug}];
            ylim([0 lim]);
            set(gca,'YDir','reverse');
            %         B = axes;
            %         set(B,'yaxislocation','right','ytick',[],'xtick',[])
            %         ylabel(B,'Right Label')
            newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Affine/MP.JN.OT/2-7/All Categories/WholeItemROIs/',presLabel{pres},'/Fixation/Raster/'];
            saveName=[newdir,'ROI_Set2-7_',titleLabel];
            
            if ~exist(newdir,'dir');
                mkdir(newdir);
            end
            
            if saveEPS
                if savePNG
                    saveas(gcf,[saveName,'.eps'],'eps')
                    title(titleLabel)
                    saveas(gcf,[saveName,'.png'],'png')
                end                
            elseif savePNG
                title(titleLabel)
                saveas(gcf,[saveName,'.png'],'png')
            else
            end
            close(gcf)
        end
    end
end

%%

% figure(101);clf(101);figure(101);hold on;
Nunits = 12;
colors = getmap(12);

% THIS WORKS
% tick_y = cell(size(trlind,1),1);
shifty = zeros(size(stim,1),1);
shift = 0;
f=permute(D{c,1},[3 2 1]);
[trial,spiketimes] = find(D{1,1}(k,:));
figure;
raster2(D{c,1},begTim:endTim,'b',shift);

for c = 1:size(D,1);
    raster(D{c,1},begTim:endTim,'b',shift);hold on;
    shifty(c,1)=shift;
    tick_y(c,1)=shift+(size(dens{c,1}{u,1},1)/2);
    shift = shift + (size(dens{c,1}{u,1},1));hold on;
end
set(gca,'ytick',tick_y{u,1}(:,1),'yticklabel',numtrl_name);
line([750 750],[shifty(1,1) (shifty(10,1)+size(dens{c,1}{u,1},1))],'LineStyle','--','Color','k','LineWidth',1);
title([filename,'__',labl{u},'__raster']);
saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Gothard Amygdala/Figures/Single Unit/Raster/',filename,'__',labl{u},'__raster','.png'], 'png')
close(gcf)


figure(101);clf(101);figure(101);hold on;
Nunits = 12;
colors = getmap(12);

% THIS WORKS
tick_y = cell(size(trlind,1),1);
for u = 1:size(trlind,1)
    shifty = zeros(size(stim,1),1);
    shift = 0;
    figure(u);
    numtrl_name=cell(size(stim,1),1);
    for c = 1:size(stim,1);
        raster(dens{c,1}{u,1},[],colors(c,:),shift);hold on;
        shifty(c,1)=shift;
        tick_y{u,1}(c,1)=shift+(size(dens{c,1}{u,1},1)/2);
        shift = shift + (size(dens{c,1}{u,1},1)+10);hold on;
        numtrl_name{c,1}=[stim_name{c,1},'_',num2str(trlnum{u,1}{c,1}),'_trls'];
    end
    set(gca,'ytick',tick_y{u,1}(:,1),'yticklabel',numtrl_name);
    line([750 750],[shifty(1,1) (shifty(10,1)+size(dens{c,1}{u,1},1))],'LineStyle','--','Color','k','LineWidth',1);
    title([filename,'__',labl{u},'__raster']);
    saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Gothard Amygdala/Figures/Single Unit/Raster/',filename,'__',labl{u},'__raster','.png'], 'png')
    clf(u)
end


%%
channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% channelcolormap(c,:)
% colorPlot=nicejet(6);
% colorPlot=hsv(6);
colorPlot=hsv(64);
% colorPlotD=brighten(-.2);

figure;
c=1;
for level=1:size(D,1)
    drug=1;
    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},colorPlot(c,:),1,smval,0,0,1,200,0)
    hold on;
    drug=2;
    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},colorPlot(c+4,:),1,smval,0,0,1,200,0)
    hold on;
    c=c+20;
end

level=3
figure;
drug=1
dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},colorPlot(c,:),1,smval,0,0,1,200,0)
hold on;
drug=2;
dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},colorPlot(c+4,:),1,smval,0,0,1,200,0)


for level=1:size(catMat,2)
    figure;
    drug=2;
    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},'b',1,smval,0,0,1,200,0)
    hold on;
    % Repeat Scenes, 2nd Presentation
    drug=1;
    dofill((begTimStr:(1/fsamp):endTimStr),D{level,drug},'r',1,smval,0,0,1,200,0)
    titleLabel=['Whole Item: ',presStr,', ',ROI_Cat_Labl{cat,2}{1,level},', ',normStr];
    title(titleLabel);
    if normArea
        if novel
            ylim([0 27])            
        else
            ylim([0 27])
        end
    else
        if novel
            ylim([0 .45])
        else
            ylim([0 .45])
        end
    end
    n1=['SL, N=',num2str(size(catMat{2,level},1))];
    n2=['OT, N=',num2str(size(catMat{1,level},1))];
    legend(n1,n2,'Location','NorthEast');
    ylabel(yStr);
    xlabel('Time From Stimulus Onset (s)');
    if ~exist(newDir,'dir');
        mkdir(newDir);
    end
    saveName=[newDir,'PSTH_Set2-7_Drug',mStr,'_',titleLabel,'.png'];
    saveas(gcf,saveName,'png')
    close(gcf)
end

%%



if psth
    for level=1:size(memMat,2)
        for drug=1:2
            figure;
            % Combine Novel Trials for Repeat & Replaced Scenes
            dofill(((1/fs):(1/fs):6),[memMat{drug,level}{1,1};memMat{drug,level}{1,1}],channelcolormap(1,:),1,smval,0,0,1,200,0)
            hold on;
            % Repeat Scenes, 2nd Presentation
            dofill(((1/fs):(1/fs):6),memMat{drug,level}{2,1},[0 1 0],1,smval,0,0,1,200,0)
            hold on;
            % Replaced Scenes, 2nd Presentation
            dofill(((1/fs):(1/fs):6),memMat{drug,level}{2,2},channelcolormap(2,:),1,smval,0,0,1,200,0)
            titleLabel=['Whole Item: ',ROI_Cat_Labl{cat,2}{level},', ',selLabl,', ',drugLabel{drug}];
            title(titleLabel);
            ylim([0 .4])
            n1=['Novel, N=',num2str(size([memMat{drug,level}{1,1};memMat{drug,level}{1,1}],1))];
            n2=['Repeat, N=',num2str(size(memMat{drug,level}{2,1},1))];
            n3=['Replaced, N=',num2str(size(memMat{drug,level}{2,2},1))];
            legend(n1,n2,n3,'Location','NorthEast');
            ylabel('Probability of Viewing ROI (%)');
            xlabel('Time From Stimulus Onset (s)');
            if ~exist(newdir,'dir');
                mkdir(newdir);
            end
            saveName=[newdir,'Memory_PSTH_Set2-7_',mStr,'_',titleLabel,'.png'];
            saveas(gcf,saveName,'png')
            close(gcf)
        end
    end
else
    for level=1:size(memMat,2)
        for drug=1:2
            D{1,drug}=sum(memMat{drug,level}{1,1},2);
            D{2,drug}=sum(memMat{drug,level}{2,1},2);
            D{3,drug}=sum(memMat{drug,level}{2,2},2);
        end
        options.stars=1;
        options.color=nicejet(size(ROI_Cat{1,1}{cat,1},2));
        options.labels={'Novel','Repeat','Replaced'};
        options.dotplot=0;
        options.test_type='ttest2';
        options.checktesttype=0;
        options.color=nicejet(2);
        options.a=.05;
        options.siglvl=1;
        options.group=1;
        barsig2(D,options)
        
        yLbl='Total Fixation Duration (ms)';
        
        titleLabel=['Whole Item: ',ROI_Cat_Labl{cat,2}{level},', ',selLabl,', ',num2str(begTimStr),'-',num2str(endTimStr),'sec'];
        title(titleLabel);

        h_legend=legend('Saline','Oxytocin','Location','NorthEast');
        set(h_legend,'FontSize',12)
        ylabel(yLbl)
        if ~exist(newdir,'dir');
            mkdir(newdir);
        end
        saveName=[newdir,'Memory_Bar_Set2-7_',mStr,'_',titleLabel,'.png'];
        saveas(gcf,saveName,'png')
        close(gcf)
        clear D
    end        
end

%%
trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,1) & include{setloop,1}(:,fileloop);
                    
                    cellMat=fix_Whole{setloop,1}{1,fileloop}(trlInd,pres);
                    xMat=zeros(size(cellMat,1),1);
                    for k=1:size(cellMat,1)
                        if ~isempty(cellMat{k,1})
                            tempMat=zeros(1,1);
                            for fixlop=1:size(cellMat{k,1},1)
                                if cellMat{k,1}(fixlop,7)<=endTim && cellMat{k,1}(fixlop,7)>=begTim
                                    tempMat(fixlop,1)=cellMat{k,1}(fixlop,9);
                                end
                            end
%                             if ~isnan(tempMat)
                                xMat(k,1)=sum(tempMat);
                                clear tempMat
%                             end
                        end
                    end
                    
                    xArea=ROI_All_area{setloop,1}(trlInd,1);
                    sceneList=ROI_Scene{setloop,1}(trlInd,1);
                    sceneListU=unique(sceneList);
                    sceneSum=nan(size(sceneListU));
                    for k=1:size(sceneListU,1);
                        ind=sceneList==sceneListU(k,1);
                        % Only calculate if there is at least one real value
                        if sum(~isnan(xMat(ind,1)))>0
                            % Only calculate if there is at least one real value
                            if sum(~isnan(xMat(ind,1)))>0
                                % Normalize looking time by ROI area if category is one
                                % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
                                if normArea
                                    if cat>2 && cat~=6
                                        sceneSum(k,1)=nansum(xMat(ind,1))./(sum(xArea(ind,:))/totalArea);
                                    else
                                        sceneSum(k,1)=nansum(xMat(ind,1));
                                    end
                                else
                                    sceneSum(k,1)=nansum(xMat(ind,1));
                                end
                            end
                        end
                    end
                    
                    clear xMat
                    if indOT{setloop,1}(1,fileloop);
                        OT=[OT;sceneSum];
                        %                         OT{1,ot}=sceneSum;
                        ot=ot+1;
                    else
                        SL=[SL;sceneSum];
                        %                         SL{1,sl}=sceneSum;
                        sl=sl+1;
                    end
                    clear xMat