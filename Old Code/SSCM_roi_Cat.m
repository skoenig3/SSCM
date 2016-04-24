

%% Calculate Looking Time in ROIs (Size,Item type,Identity,Eyes
% Visible,Age,Sex)
vardir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';

load([vardir 'ROI_Cat_Set1_xy_labl']);
load([vardir 'SSCM_Set1_MP.IW.JN_eyedata.mat']);
load([vardir 'Eyedatpixlr_S1_MP.JN.IW.mat']);
load([vardir 'SSCM90_Set1_cnd.mat']);
% load([savdir 'ROI_Cat_ID.mat']);
clear eyedat

% 
% ROI_IDs=ROI_ID;
% ROI_ID_Labl=ROI_Cat_Labl(1:2,1:2);
% ROI_ID_Labl(3:5,1:2)=ROI_Cat_Labl(4:6,1:2);
% ROI_Cat_Labl=ROI_ID_Labl;
% ROI_M_IDs=ROI_ID;
% 
% 
% %save([savdir 'ROI_Cat_ID.mat'],'ROI_Cat_xy','ROI_Cat_Labl','ROI_M_IDs','ROI_M_IDs_Labl');
% 
% ROI_ID=ROI_Cat{3,1};
% ROI_Cat2=cell(5,1);
% ROI_Cat2{1}=ROI_Cat{1};
% ROI_Cat2{2}=ROI_Cat{2};
% ROI_Cat2(3:5,1)=ROI_Cat(4:6,1);
% ROI_Cat=ROI_Cat2;
% clear ROI_Cat2
% 
% ROI_Cat_xy2=cell(5,1);
% ROI_Cat_xy2{1}=ROI_Cat_xy{1};
% ROI_Cat_xy2{2}=ROI_Cat_xy{2};
% ROI_Cat_xy2(3:5,1)=ROI_Cat_xy(4:6,1);
% ROI_Cat_xy=ROI_Cat_xy2;
% clear ROI_Cat_xy2



%Novel Presentation
tic
roi_Cat=cell(size(ROI_Cat));
roi_Cat_eyedat=cell(size(ROI_Cat));
for catloop=1:size(ROI_Cat,1);
    for level=1:size(ROI_Cat{catloop},2);
        for fileloop=1:length(M);
            for roi=1:length(ROI_Cat{catloop}{level});
                scene=str2double(ROI_Cat{catloop}{level}{roi,1}(13:14));
                % Novel Presentation (+1010, +1011 for 2nd Presentation)
                cond=(((scene-1)*2)+1010);
                ind=find(cnd==cond);
                if ~isempty(ind)
                    if size(eyedatpixlr{ind,fileloop},2)>=2000
                        % Select ROI
                        regions=ROI_Cat_xy{catloop}{level}{roi,1};
                        % Create logical index of timepoints when eye
                        % position coincided with ROI
                        roi_ind=ismember(eyedatpixlr{ind,fileloop}',regions,'rows');
                        roi_Cat{catloop,1}{1,level}{1,fileloop}(roi,:)=roi_ind(1:2000,1)';
                        % Select eye positions that coincided with ROI
                        x=eyedatpixlr{ind,fileloop}(1,:);
                        y=eyedatpixlr{ind,fileloop}(2,:);
                        roi_eyedat=[];
                        roi_eyedat(1,:)=x(roi_ind');
                        roi_eyedat(2,:)=y(roi_ind');
                        roi_Cat_eyedat{catloop}{level}{roi,fileloop}=roi_eyedat;
                    end
                end
            end
        end
    end
end
    toc
    
%%
vardir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
%save([vardir 'roi_Cat+Labl.mat'],'roi_Cat','ROI_Cat_Labl','M','Set','cnd','cndFile','itmFile','itmfil');

%%
vardir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
load([vardir 'roi_Cat+Labl.mat']);
ROI_Cat(:,2:3)=ROI_Cat_Labl;

% Plot Prob of Viewing ROI Averaged Across All Trials for All Monkeys in
% Novel Scenes for Each Level of Each Category 
% (Size,Item,Type,Eyes,Age,Sex)
colors = [0 0 1;1 0 0;.1 .9 .1;1 0 1];
for catloop=1:length(roi_Cat);
    figure;
    label=[];
    for level=1:length(roi_Cat{catloop,1});
        if strfind(ROI_Cat_Labl{catloop,2}{1,level},'10000')
          dofill((.005:.005:10),([+roi_Cat{catloop,1}{1,level}{1,1};+roi_Cat{catloop,1}{1,level}{1,2};+roi_Cat{catloop,1}{1,level}{1,3}]/(10000/(600*800))), colors(level,:),1,60,0,0,1,200,1)             
        elseif strfind(ROI_Cat_Labl{catloop,2}{1,level},'5000')
          dofill((.005:.005:10),([+roi_Cat{catloop,1}{1,level}{1,1};+roi_Cat{catloop,1}{1,level}{1,2};+roi_Cat{catloop,1}{1,level}{1,3}]/(5000/(600*800))), colors(level,:),1,60,0,0,1,200,1)             
        elseif strfind(ROI_Cat_Labl{catloop,2}{1,level},'2000')    
          dofill((.005:.005:10),([+roi_Cat{catloop,1}{1,level}{1,1};+roi_Cat{catloop,1}{1,level}{1,2};+roi_Cat{catloop,1}{1,level}{1,3}]/(2000/(600*800))), colors(level,:),1,60,0,0,1,200,1)             
        else
          dofill((.005:.005:10),[+roi_Cat{catloop,1}{1,level}{1,1};+roi_Cat{catloop,1}{1,level}{1,2};+roi_Cat{catloop,1}{1,level}{1,3}], colors(level,:),1,60,0,0,1,200,1)
%         plot((.001:.005:10),smooth((100*nanmean([+roi_Cat{catloop,1}{1,level}{1,1};+roi_Cat{catloop,1}{1,level}{1,2};+roi_Cat{catloop,1}{1,level}{1,3}])),20),'LineStyle' , '-' , 'Color' , colors(level,:),'LineWidth',1);
        end        
        hold on;
        if level~=length(roi_Cat{catloop,1});
            label=[label ROI_Cat_Labl{catloop,2}{1,level} ', '];
        else
            label=[label ROI_Cat_Labl{catloop,2}{1,level}];
        end
        legendVar{level}=[ROI_Cat_Labl{catloop,2}{1,level},'   n=',num2str(size(roi_Cat{catloop,1}{1,level}{1,1},1)*length(roi_Cat{catloop,1}{1,level}))];
    end
    legend(legendVar,'Location','Best')
    xlabel('Time From Stimulus Onset (s)')
%     ylabel('Probability ofViewing ROI (%)')
    ylabel('ROI Hit Rate (Hz)')
    title([ROI_Cat_Labl{catloop,1}(1:end) ': ' label])
%     newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/'];
%     if ~exist(newdir,'dir');
%         mkdir(newdir);
%     end
%     %saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/','ROI_Set',num2str(Set),'_',ROI_Cat_Labl{catloop,1},'_MP.IW.JN','_ProbView','.png'], 'png')
%     %saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/ROI_Set',num2str(Set),'_',ROI_Cat_Labl{catloop,1},'_MP.IW.JN','_ProbView','.eps'], 'eps')
%     close all
end

.07/(10000/(600*800))
.04/(5000/(600*800))
.02/(2000/(600*800))

%% Plot Prob of Viewing ROI Averaged Across All Trials for Each Monkey in
% Novel Scenes for Each Level of Each Category 
% Separate Figures for Each Monkey
% (Size,Item,Type,Eyes,Age,Sex)

colors = [0 0 1;1 0 0;.1 .9 .1;1 0 1];
for catloop=1:length(roi_Cat);
%     maxM=zeros(length(M),length(roi_Cat{catloop,1}));
%     for monkloop=1:length(M);
%         for level=1:length(roi_Cat{catloop,1});
%             maxM(monkloop,level)=max(smooth((100*nanmean(+roi_Cat{catloop,1}{1,level}{1,monkloop})),20));
%         end
%     end
%     maxM=round(1.1*max(max(maxM)));
    for monkloop=1:length(M);
        figure;
        label=[];
        maxMum=zeros(length(roi_Cat{catloop,1}),1);
        for level=1:length(roi_Cat{catloop,1});
            dofill((.005:.005:10),+roi_Cat{catloop,1}{1,level}{1,monkloop}, colors(level,:),1,60,0,0,1,200,1)
            maxMum(level,1)=ans(2);
            %Plot probability of a hit being recorded on an roi,
            %sampled at 200 Hz and smoothed in bins of 100 ms
%             plot((.005:.005:10),smooth((100*nanmean(+roi_Cat{catloop,1}{1,level}{1,monkloop})),20),'LineStyle' , '-' , 'Color' , colors(level,:),'LineWidth',1);
%             dofill(.005:.005:10,+roi_Cat{catloop,1}{1,level}{1,monkloop},1,1,0,0,0);
            hold on;
            if level~=length(roi_Cat{catloop,1});
                label=[label ROI_Cat_Labl{catloop,2}{1,level} ', '];
            else
                label=[label ROI_Cat_Labl{catloop,2}{1,level}];
            end
            legendVar{level}=[ROI_Cat_Labl{catloop,2}{1,level},'   n=',num2str(size(roi_Cat{catloop,1}{1,level}{1,monkloop},1))];
        end 
        maxMum=max(maxMum);
        legend(legendVar,'Location','Best')
        ylim([0 maxMum]);
        xlabel('Time From Stimulus Onset (s)')
        ylabel('Probability of Viewing ROI (%)')
        title([M{1,monkloop},', ' ROI_Cat_Labl{catloop,1}(1:end),': ',label])
        newdir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/'];
        if ~exist(newdir,'dir');
            mkdir(newdir);
        end
        %saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/','ROI_Set',num2str(Set),'_',ROI_Cat_Labl{catloop,1},'_',M{1,monkloop},'_ProbView','.png'], 'png')
        %saveas(gcf,['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Figures/Set',num2str(Set),'/',ROI_Cat_Labl{catloop,1},'/ROI_Set',num2str(Set),'_',ROI_Cat_Labl{catloop,1},'_',M{1,monkloop},'_ProbView','.eps'], 'eps')
        close all
    end
end

 figure;dofill((.001:.005:10),y,'red',1,60)
 sampling_rate = 200;
 spikes = rand(100,2000);
 time = 0.001:0.001:2;
 figure;
 dofill((.001:.005:10),y,'red',1,60,0,0,1,sampling_rate,1)

% colors = {[0 0 1],[1 0 0],[.1 .9 .1];[0 0 .8],[.8 0 0],[.1 .7 .1];[0 0 .6],[.6 0 0],[.1 .5 .1]};

%% Area - Item Type
vardir='R:\Buffalo Lab\eblab\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
load([vardir 'roi_Cat+Labl.mat']);
% roi_Cat(:,2:3)=ROI_Cat_Labl;
load([vardir 'ROI_Cat_Set1.mat']);

ROI_Size_Item=cell(length(ROI_Cat{1,1}),length(ROI_Cat{2,1}));
for sizeloop=1:length(ROI_Cat{1,1});
    for itemloop=1:length(ROI_Cat{2,1});
        ROI_Size_Item{sizeloop,itemloop}=intersect(ROI_Cat{1,1}{1,sizeloop},ROI_Cat{2,1}{1,itemloop});
    end
end

