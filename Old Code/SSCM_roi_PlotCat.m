% %% Plot Bar Graphs of Looking Time at Category: Whole-Item ROIS
% channelcolormap = [0.75 0 0;0 0 1;0 1 0;0.44 0.19 0.63;0 0.13 0.38;0.5 0.5 0.5;1 0.75 0;1 0 0;0.89 0.42 0.04;0.85 0.59 0.58;0.57 0.82 0.31;0 0.69 0.94;1 0 0.4;0 0.69 0.31;0 0.44 0.75];
% vardir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
% inddir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
% load([inddir 'SSCM_ROI_Include.mat']);
% load([inddir 'SSCM_indOT.mat']);
% load([inddir 'SSCM_roi_All_ReCal_130112_00h40m.mat']);
% load([inddir 'SSCM_Set1.2.3.4.5.6.7_ROI_All_xy.mat']);
% load([inddir 'SSCM_Set1.2.3.4.5.6.7_indTrl_ROI_Cat.mat']);
% clear ROI_All_xy

%% Plot Probability of Viewing Monkeys & Objects on OT & SL
% Novel Only
smval=40;
%Sampling Frequency
fsamp=200;
%Total Area of Image in Pixels
totalArea=800*600;

% Normalize Viewing Time to ROI Area in Unbalanced Conditions?
normArea=true;

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


% What Time Period??????
begTimStr=0;
if begTimStr==0
    begTim=1;
else
    begTim=begTimStr*fsamp;
end
endTimStr=10;
endTim=endTimStr*fsamp;

cat=2;
clear catMat
clear trlInd
for d=1:2
    for k=1:size(ROI_Cat_Labl{cat,2},2)
        start(d,k)=1;
    end
end
for level=1:size(ROI_Cat_Labl{cat,2},2);
    for setloop=2:7
        for fileloop=1:length(roi_All_1{setloop,1})
            if ~isempty(strfind(mStr,Sets{setloop,1}{1,fileloop}(1:2)))
                
                trlInd=ROI_Cat{setloop,1}{cat,1}(:,level) & indTrl{setloop,1}(:,1) & include{setloop,1}(:,fileloop);
                xMat=roi_All_1{setloop,1}{1,fileloop}(trlInd,begTim:endTim);
                
                xArea=ROI_All_area{setloop,1}(trlInd,1);                                    
                sceneList=ROI_Scene{setloop,1}(trlInd,1);
                sceneListU=unique(sceneList);
                sceneSum=zeros(size(sceneListU,1),size((begTim:endTim),2));
                for k=1:size(sceneListU,1);
                    ind=sceneList==sceneListU(k,1);
                    % Normalize looking time by ROI area if category is one
                    % that isn't balanced across scenes (Area of ROI,Item Type or Gaze)
%                     if cat>2 && cat~=6                        
                        if normArea
                            sceneSum(k,begTim:endTim)=((mean(xMat(ind,begTim:endTim),1)).\(sum(xArea(ind,:))\totalArea));
                        else
                            sceneSum(k,begTim:endTim)=mean(xMat(ind,begTim:endTim),1);
                        end
%                     else
%                         sceneSum(k,begTim:endTim)=mean(xMat(ind,begTim:endTim),1);                        
%                     end
                end
                
                if indOT{setloop,1}(1,fileloop);
                    drug=1;
                    catMat{drug,level}((start(drug,level):(start(drug,level)+(size(sceneSum,1)-1))),:)=sceneSum;
                    start(drug,level)=start(drug,level)+size(sceneSum,1);
                else
                    drug=2;
                    catMat{drug,level}((start(drug,level):(start(drug,level)+(size(sceneSum,1)-1))),:)=sceneSum;
                    start(drug,level)=start(drug,level)+size(sceneSum,1);
                end
                clear sceneSum
            else
            end
        end
    end
end

%% Plot OT & SL For Each Level
drugLabel={'OT','SL'};
for level=1:size(catMat,2)
    figure;
    % Combine Novel Trials for Repeat & Replaced Scenes
    dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{2,level},'b',1,smval,0,0,1,200,0)
    hold on;
    % Repeat Scenes, 2nd Presentation
    dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{1,level},'r',1,smval,0,0,1,200,0)
    titleLabel=['Whole Item: ',drugLabel{drug},' ',ROI_Cat_Labl{cat,2}{1,level}];
    title(titleLabel);
%     ylim([0 .4])
%     n1=['Novel, N=',num2str(size([memMat{drug,level}{1,1};memMat{drug,level}{1,1}],1))];
%     n2=['Repeat, N=',num2str(size(memMat{drug,level}{2,1},1))];
%     n3=['Replaced, N=',num2str(size(memMat{drug,level}{2,2},1))];
%     legend(n1,n2,n3,'Location','NorthEast');
    ylabel('Probability of Viewing ROI (%)');
    xlabel('Time From Stimulus Onset (s)');
%     if ~exist(newdir,'dir');
%         mkdir(newdir);
%     end
%     %saveName=[newdir,'Memory_PSTH_Set2-7_',mStr,'_',titleLabel,'.png'];
%     %saveas(gcf,%saveName,'png')
%     close(gcf)
end

%% Plot OT & SL Separately
drugLabel={'OT','SL'};
for drug=1:2
    figure;
    % Combine Novel Trials for Repeat & Replaced Scenes
    dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,1},'b',1,smval,0,0,1,200,0)
    hold on;
    % Repeat Scenes, 2nd Presentation
    if size(catMat,2)==3
        dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,2},'r',1,smval,0,0,1,200,0)
        hold on;
        dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,3},'g',1,smval,0,0,1,200,0)
    else
        dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,2},'r',1,smval,0,0,1,200,0)
    end
     
       
        
    titleLabel=['Whole Item: ',drugLabel{drug},' ',ROI_Cat_Labl{cat,1}];
    title(titleLabel);
%     ylim([0 .4])
%     n1=['Novel, N=',num2str(size([memMat{drug,level}{1,1};memMat{drug,level}{1,1}],1))];
%     n2=['Repeat, N=',num2str(size(memMat{drug,level}{2,1},1))];
%     n3=['Replaced, N=',num2str(size(memMat{drug,level}{2,2},1))];
%     legend(n1,n2,n3,'Location','NorthEast');
    ylabel('Probability of Viewing ROI (%)');
    xlabel('Time From Stimulus Onset (s)');
%     if ~exist(newdir,'dir');
%         mkdir(newdir);
%     end
%     %saveName=[newdir,'Memory_PSTH_Set2-7_',mStr,'_',titleLabel,'.png'];
%     %saveas(gcf,%saveName,'png')
%     close(gcf)
end

%% Plot OT & SL Together
newDir='\Users\DrewSolyst\Documents\Buffalo Rotation\Scene Manipulation\Figures\Affine\MP.JN.OT\2-7\All Categories\WholeItemROIs\Attention\PSTH\';

fs=200;
frate=0;
smval=40;
figure;
drug=2;
% Objects SL
dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,2},'b',1,smval,0,0,1,fs,frate)
hold on;
% Objects OT
drug=1;
dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,2},'r',1,smval,0,0,1,fs,frate)
hold on;
drug=2;
% Monkeys SL
dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,1},'b',1,smval,0,0,1,fs,frate)
hold on;
drug=1;
% Monkeys OT
dofill(((1\fsamp):(1\fsamp):endTimStr),catMat{drug,1},'r',1,smval,0,0,1,fs,frate)
titleLabel='Probability of Viewing Monkeys & Objects in Novel Scenes, OT & SL';
title(titleLabel);
ylim([0 1])
n1=['Saline, N=',num2str(size(catMat{1,1},1))];
n2=['Oxytocin, N=',num2str(size(catMat{2,1},1))];
legend(n1,n2,'Location','NorthEast');
ylabel('Probability of Viewing ROI (%)');
xlabel('Time From Stimulus Onset (s)');
% if ~exist(newDir,'dir');
%     mkdir(newDir);
% end
% %saveName=[newDir,'PSTH_Set2-7_',mStr,'_',titleLabel,'.png'];
% %saveas(gcf,%saveName,'png')
% close(gcf)

