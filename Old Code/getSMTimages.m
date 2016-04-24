cndSceneList=cell(1,size(trlInd,2));
cndTrialList=cell(1,size(trlInd,2));
for k=1:size(trlInd,2)
    cndTrialList{k}=find(trlInd(:,k)==1);
    cndSceneList{k}=round((find(trlInd(:,k)==1))/2);
end

% Include ROI_M,ROI_M_Di,ROI_M_Av,ROI_O in RepT
% ROI_cndInd=zeros(length(ROI_Set1),size(trlInd,2));
% ROI_cndInd(:,1)=[1;1;1;0;0;0;1;0;0;0];
% ROI_cndInd(:,2)=[1;1;1;0;0;1;1;0;0;1];
% ROI_cndInd(:,3)=[1;1;1;0;0;1;1;0;0;1];
% ROI_cndInd(:,4)=[0;0;0;0;0;0;0;1;0;0];
% ROI_cndInd(:,5)=[0;0;0;0;0;0;0;0;1;0];
% ROI_cndInd(:,6)=[0;0;0;1;0;0;0;0;0;0];
% ROI_cndInd(:,7)=[0;0;0;0;1;0;0;0;0;0];
% ROI_cndInd=logical(ROI_cndInd);

% No Overlap of ROIs
ROI_cndInd=zeros(length(ROI_Set1),size(trlInd,2));
ROI_cndInd(:,1)=[1;1;1;0;0;0;1;0;0;0];%Novel 
ROI_cndInd(:,2)=[0;0;0;0;0;1;0;0;0;1];%Repeat: 1st Presentation
ROI_cndInd(:,3)=[1;1;1;0;0;1;1;0;0;1];%Repeat: 2nd Presentation
ROI_cndInd(:,4)=[0;0;0;0;0;0;0;1;0;0];%Object Replaced: 1st Presentation
ROI_cndInd(:,5)=[0;0;0;0;0;0;0;0;1;0];%Object Replaced: 2nd Presentation
ROI_cndInd(:,6)=[0;0;0;1;0;0;0;0;0;0];%Monkey Replaced: 1st Presentation
ROI_cndInd(:,7)=[0;0;0;0;1;0;0;0;0;0];%Monkey Replaced: 2nd Presentation
ROI_cndInd=logical(ROI_cndInd);

ROI_cndName=cell(1,size(ROI_cndInd,2));
for k=1:size(ROI_cndInd,2);
    ROI_cndName{1,k}=ROI_Set1Name(ROI_cndInd(:,k))';
end

ROI_Set1Dir={'/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M_Av/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M_Di/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M_RepL_CR1/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M_RepL_CR2/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_M_RepT/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_O/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_O_RepL_CR1/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_O_RepL_CR2/';...
    '/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set1/ROIs/FINAL ROI/ROI_O_RepT/'};

save('/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Processed Data/EyeData_Set_1_Monkeys-4_ROI.mat','M','ROI_Set1','ROI_Set1Dir','ROI_Set1Name','ROI_cndInd','S','cndFileList','cndList','cndSceneList','cndFileList','eyedat','itmFileList','trlInd','trlIndLabl','xcal','ycal','files');

%% TEST
subjloop=1;
cndloop=1;
trlloop=6;
roiloop=2;
cndROI=ROI_Set1(ROI_cndInd(:,cndloop),S(setloop));
cndROIdir=ROI_Set1Dir(ROI_cndInd(:,cndloop),S(setloop));
            

if cndSceneList{1,cndloop}(trlloop,1)<10
    scene=['0' num2str(cndSceneList{1,cndloop}(trlloop,1))];
else
    scene=num2str(cndSceneList{1,cndloop}(trlloop,1));
end
ind=strfind(cndROI{roiloop,1},scene);
isFound=~cellfun('isempty',ind);
image= cndROI{roiloop,1}(isFound);
I = imread([cndROIdir{roiloop} image{1:end}]);
I=rgb2gray(I);
[Y,X]=find(I > 0);
regions=[X 600-Y];
regions=regions';

figure;scatter(round((eyedat{1,1}{11,1}(1,:)*24)+400),round((eyedat{1,1}{11,1}(2,:)*24)+300));axis([0 800 0 600])
figure;scatter(regions(1,:),regions(2,:));axis([0 800 0 600])
figure;scatter(trialspix{trlloop,1}(1,:),trialspix{trlloop,1}(2,:));axis([0 800 0 600])

for i=1:size(trials{trlloop,1},2);
    inRegion_hor=find(regions(1,:)==trialspix{trlloop,1}(1,i));
    inRegion_vrt=find(regions(2,:)==trialspix{trlloop,1}(2,i));
    inRegion=intersect(inRegion_hor,inRegion_vrt);
    if ~isempty(inRegion)
        viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(1,i)=1;
    else
    end
end

looktime{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=sum(viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop})*5;

%%
load('/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Processed Data/EyeData_Set_1_Monkeys-4_ROI.mat');

%%

% construct filter for eye velocity measure
fltord = 60;
lowpasfrq = 25;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);


transitionsIn=cell(length(S),size(files,2));
transitionsOut=cell(length(S),size(files,2));
looktime=cell(length(S),size(files,2));
fixationsIn=cell(length(S),size(files,2));
fixationsOut=cell(length(S),size(files,2));
fixationsFirst=cell(length(S),size(files,2));
fixationsLast=cell(length(S),size(files,2));
% Allocate cells of nSets (1) x nSubjects (4)
viewROI=cell(length(S),size(files,2));

tic
for setloop=1:length(S);
    for subjloop=1:size(files,2);
        transitionsIn{setloop,subjloop}=cell(1,size(trlInd,2));
        transitionsOut{setloop,subjloop}=cell(1,size(trlInd,2));
        looktime{setloop,subjloop}=cell(1,size(trlInd,2));
        fixationsIn{setloop,subjloop}=cell(1,size(trlInd,2));
        fixationsOut{setloop,subjloop}=cell(1,size(trlInd,2));
        fixationsFirst{setloop,subjloop}=cell(1,size(trlInd,2));
        fixationsLast{setloop,subjloop}=cell(1,size(trlInd,2));
        % For each subject, Allocate cells of nConditions (7)
        viewROI{setloop,subjloop}=cell(1,size(trlInd,2));
        for cndloop=1:size(trlInd,2);
            % Select condition to analyze (Novel,Repeat 1st,etc)
            trials=eyedat{setloop,1}(trlInd(:,cndloop),subjloop);
            trialspix=cell(size(trials));
            for k=1:length(trials);
                trialspix{k}=zeros(2,size(trials{k},2));
            end
            
            % trialspixvrt had +400 instead of +300!!! Mundane details
            % Michael!
            for k=1:size(trials,1);
                for i=1:size(trials{k,1},2);
                    trialspix{k,1}(1,i)=round(((trials{k,1}(1,i))*24)+400);
                    trialspix{k,1}(2,i)=round(((trials{k,1}(2,i))*24)+300);
                end
            end
                    
            figure;scatter(trialspix{trlloop,1}(1,:),trialspix{trlloop,1}(2,:));axis([0 800 0 600])
            figure;scatter(trials{trlloop,1}(1,:),trials{trlloop,1}(2,:));
            axis([-33.333 33.333 -15 15])
            
            cndROI=ROI_Set1(ROI_cndInd(:,cndloop),S(setloop));
            cndROIdir=ROI_Set1Dir(ROI_cndInd(:,cndloop),S(setloop));
            
            transitionsIn{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            transitionsOut{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            looktime{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            fixationsIn{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            fixationsOut{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            fixationsFirst{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            fixationsLast{setloop,subjloop}{1,cndloop}=zeros(length(trials),length(cndROI));
            
            % For subject of subjloop and condition of cndloop,
            % Allocate cells of nTrials by nRois (cndROI)
            viewROI{setloop,subjloop}{1,cndloop}=cell(length(trials),length(cndROI));

            for trlloop=1:length(trials);
                % differentiate and multiply with sampling rate to get velocity as deg/sec
                resampx=resample(trials{trlloop}(1,:),5,1);
                resampy=resample(trials{trlloop}(2,:),5,1);
                
                % calculate velocity
                x_v= diff(filtfilt(flt,1, resampx)) .* 1000;
                y_v= diff(filtfilt(flt,1, resampy)) .* 1000;
                
                % combine x- and y-velocity to get overall eye velocity
                vel = abs(complex(x_v,y_v));
                
                % lim = threshold for detecting saccade
                lim = 40;
                sacbeg = find(diff(vel > lim) > 0);
                sacend = find(diff(vel > lim) < 0);
                
                if vel(end)>=lim
                    if vel(1)<lim
                        tempbeg=[1 sacend];
                        tempend=sacbeg;
                    else
                        tempbeg=sacend;
                        tempend=sacbeg;
                    end
                else
                    if vel(1)<lim
                        tempbeg=[1 sacend];
                        tempend=[sacbeg length(vel)];
                    else
                        tempbeg=sacend;
                        tempend=[sacbeg length(vel)];
                    end
                end
                
                % Load ROIs and select appropriate ROI
                for roiloop=1:length(cndROI);
                    % Each roiloop has a different # of trials !!!!!!!
                    % For subject of subjloop, condition of cndloop, and
                    % ROI of roiloop, Allocate zeros of n
                    viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}=zeros(1,size(trials{trlloop,1},2));

                    if cndSceneList{1,cndloop}(trlloop,1)<10
                        scene=['0' num2str(cndSceneList{1,cndloop}(trlloop,1))];
                    else
                        scene=num2str(cndSceneList{1,cndloop}(trlloop,1));
                    end
                    ind=strfind(cndROI{roiloop,1},scene);
                    isFound=~cellfun('isempty',ind);
                    image= cndROI{roiloop,1}(isFound);
                    I = imread([cndROIdir{roiloop} image{1:end}]);
                    I=rgb2gray(I);
                    [Y,X]=find(I > 0);
                    regions=[X 600-Y];
                    regions=regions';
                    clear Y,clear X,clear I,clear isFound
                    
                   
                    % Enter a "1" at every timepoint when eye position (in pixels)
                    % falls within coordinates specified in the loaded ROI
                    % REALLY REALLY REALLY SLOW!!!
%                     for i=1:size(trials{trlloop,1},2);
%                         for k=1:length(regions);
%                             if regions(1,k)==trialspix{trlloop,1}(1,i) && regions(2,k)==trialspix{trlloop,1}(2,i)
%                                 viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(i)=1;
%                             else
%                             end
%                         end
%                     end
                    
                    %Enter a "1" at every timepoint when eye position (in pixels)
                    %falls within coordinates specified in the loaded ROI
                    % ----To Plot the parts of the ROI that the Monkey Viewed---
                    % store: regions(inRegion)
                    for i=1:size(trials{trlloop,1},2); 
                        inRegion_hor=find(regions(1,:)==trialspix{trlloop,1}(1,i));
                        inRegion_vrt=find(regions(2,:)==trialspix{trlloop,1}(2,i));
                        inRegion=intersect(inRegion_hor,inRegion_vrt);
                        if ~isempty(inRegion)
                            viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(1,i)=1;
                        else
                        end
                    end
                    
                    %Calculate time spent viewing ROI as total number of
                    %instances viewing CR multiplied by 5 (1 sample every
                    % 5ms) to get looking time in ms
                    looktime{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=sum(viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop})*5;
                    
                    % calculate fixations in CR & transitions into/out of CR
                    fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=0;% # of fixations inside ROI
                    fixationsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=0; % # of total fixations for trial
                    transitionsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=0; % # of transitions into ROI
                    transitionsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=0; % # of transitions out of ROI
                    
                    for fixlop=1:length(tempbeg)
                        fixhor=round((mean(resampx(tempbeg(fixlop):tempend(fixlop)))*24)+400); % fixation position in pixels, horizontal
                        fixvrt=round((mean(resampy(tempbeg(fixlop):tempend(fixlop)))*24)+300); % fixation position in pixels, vertical
                        if fixlop~=1
                            fixhor_lasttrl=round(mean((resampx(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+400); % previous fixation position, horizontal
                            fixvrt_lasttrl=round(mean((resampy(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+300); % previous fixation position, horizontal
                        end
                        fixout=1;
                        
                        % Find fixations and transitions in CR
                        if inROI(fixhor,fixvrt,regions);
                            fixout=0;
                            if fixlop==1
                                fixationsFirst{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=1;
                                fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1;
                            elseif fixlop==length(tempbeg);
                                fixationsLast{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=1;
                                fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1;
                            else
                                fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=fixationsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1; %add 1 to # of fixations inside area of interest
                            end
                            if fixlop~=1
                                if ~inROI(fixhor_lasttrl,fixvrt_lasttrl,regions);
                                    transitionsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=transitionsIn{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1;
                                end
                            end
                        end
                        if fixout==1
                            if fixlop~=1
                                if inROI(fixhor_lasttrl,fixvrt_lasttrl,regions);
                                    transitionsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=transitionsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1;
                                end
                            end
                        end
                        if fixhor>=0 && fixhor<=800
                            if fixvrt>=0 && fixvrt<=600
                                fixationsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)=fixationsOut{setloop,subjloop}{1,cndloop}(trlloop,roiloop)+1;
                            end
                        end
                    end
                end
            end
        end
    end
end
toc
%% Remove points outside of screen

%---remove values outside of 800x600 pixel screen---%
for setloop=1:length(S);
    
    
    for subjloop=1:size(eyedat{1,1},2);
        
        for cndloop=1:size(trlInd,2);
            trials=eyedat{setloop,1}(trlInd(:,cndloop),subjloop);
            cndROI=ROI_Set1(ROI_cndInd(:,cndloop),S(setloop));

            for trlloop=1:length(trials);
                
                for roiloop=1:length(cndROI);
                    
                    
                    x = (trials{trlloop,1}(1,:)*24)+400;
                    y = (trials{trlloop,1}(2,:)*24)+300;
                    badx = find(x < 1 | x > 800);
                    x(badx) = []; y(badx) = [];
                    bady = find(y < 1 | y > 600);
                    x(bady) = []; y(bady) = [];
                    viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(badx) =[];
                    viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(bady) =[];
                end
                
            end
            
        end
    end
    
end

%% Error
% The first and second presentation of scene 41 is switched in viewROI
% Can't swap viewROI data because ROI's are different
temp=cell(1,4);
for k=1:4
    temp{1,k}=viewROI{1,2}{1,1}{41,k};
end
for k=1:4
    viewROI{1,2}{1,1}{41,k}=viewROI{1,2}{1,4}{13,k};
end
for k=1:4
    viewROI{1,2}{1,1}{13,k}=temp{1,k};
end

% Replace erroneous trials in viewROI with zeros
for k=1:4
    viewROI{1,2}{1,1}{41,k}=zeros(1,2000);
    viewROI{1,3}{1,1}{41,k}=zeros(1,2000);
    viewROI{1,2}{1,3}{13,k}=zeros(1,1200);
    viewROI{1,3}{1,3}{13,k}=zeros(1,1200);    
end



%% Analyze and Plot the Data: Average all (including Timmy)
cndROI=ROI_Set1(ROI_cndInd(:,cndloop),S(setloop));

roi_M=zeros(90,4);
cndloop=1;
roiloop=3;

% plotMat=zeros(length(viewROI{setloop,subjloop}{1,cndloop})*size(viewROI,2),plotTime);
for setloop=1:length(S);
%     for subjloop=1:size(viewROI,2);
    for subjloop=2;
        if size(viewROI{setloop,subjloop}{1,cndloop}{1,1},2)<1970
            % viewROI{1,3}{1,4}(17&20,1) only go to 1970 samples            
            plotTime=1173;
        else
            % viewROI{1,2}{1,7}(20&27,1) only go to 1197 samples
            % viewROI{1,3}{1,7}(7,1) only go to 1173 samples
            plotTime=1970;
        end
        plotMat=zeros(length(viewROI{setloop,subjloop}{1,cndloop})*size(viewROI,2),plotTime);
        offset=0;
        for trlloop=1:length(viewROI{setloop,subjloop}{1,cndloop});
            % The first and second presentation of scene 41 is switched in
            % viewROI for Subj 2 and 3
            plotMat(trlloop+offset,:)=viewROI{setloop,subjloop}{1,cndloop}{trlloop,roiloop}(1:plotTime);
        end
        offset=offset+length(viewROI{setloop,subjloop}{1,cndloop});
    end
end
% Remove Trials 41 of Subj 2 and 3 from 
plotMat((41+90),:)=[];
plotMat((41+180),:)=[];


%Plot Probability of Viewing Monkeys or Object in Novel Scenes
% figure;
% plot(smooth(plotMat,100));
figure;
plot(density(plotMat,100,'gauss'));
axis([0 2000 0 40])


figure;
plot((.001:.005:9.85),density(plotMat,100,'gauss'));
axis([.001 9.85 0 .75])

axis([.001 9.85 0 .75])

plot((.001:.005:plotTime),density(plotMat,100,'gauss'),'b');
axis([.001 10 0 .75])

legend(['Monkeys   ','n = ' num2str(size(roi_Mp,1));'Objects   ','n = ', num2str(size(roi_Op,1))],'Location','Best')
xlabel('Time (s)')
ylabel('Probability of Viewing ROI')
title('Probability of Viewing Monkeys or Objects in Novel Scenes')


%%
plotMat=zeros(length(viewROI{setloop,subjloop}{1,cndloop})*size(viewROI,2),plotTime);

%%
load('SSCM_ROI_viewROI_plotMat.mat');




%%
load('SCM.mat');
%% TEST
fixhortest=zeros(1,length(tempbeg));
fixvrttest=zeros(1,length(tempbeg));
for fixlop=1:length(tempbeg)
    fixhortest(fixlop)=(mean(resampx(tempbeg(fixlop):tempend(fixlop)))*24)+400; % fixation position in pixels, horizontal
    fixvrttest(fixlop)=(mean(resampy(tempbeg(fixlop):tempend(fixlop)))*24)+300; % fixation position in pixels, vertical
    
end
fixtest=round([fixhortest;fixvrttest]);
figure;scatter(fixtest(1,:),fixtest(2,:));
axis([0 800 0 600]);


eyedatpixl=zeros(2,length(trials{30,1}));
for k=1:length(trials{30,1});
    eyedatpixl(1,k)=(trials{30,1}(1,k)*24)+400;
    eyedatpixl(2,k)=(trials{30,1}(2,k)*24)+300;
end


figure;scatter(eyedatpixl(1,:),eyedatpixl(2,:));
axis([0 800 0 600]);

figure;plot(eyedatpixl(1,:),eyedatpixl(2,:));
axis([0 800 0 600]);


figure;scatter(fixhortest,fixvrttest);
axis([0 800 0 600]);

figure;scatter(regions(1,:),regions(2,:));
axis([0 800 0 600]);

figure;scatter(resampx,resampy);
axis([0 800 0 600]);

figure;scatter(trials{30,1}(1,:),trials{30,1}(2,:));
axis([0 800 0 600]);


figure;plot(vel);hold on;
line([0 7000],[40 40],'LineStyle' , '--' , 'Color' , 'k','LineWidth' , 1);

% test=intersect([fixhortest;fixvrttest],regions,'rows');
% fixtest=round([fixhortest;fixvrttest]);
% fixIn=0;
% fixInInd=[];
% fixOutInd=[];
% for k=1:length(fixtest);
%    if ismember(fixtest(:,k),regions)
%        fixIn=fixIn+1;
%        fixInInd=[fixInInd,k];
%    else
%        fixOutInd=[fixOutInd,k];
%    end
% end
fixInPlot=fixtest(:,fixInInd);
fixOutPlot=fixtest(:,fixOutInd);
figure;scatter(fixtest(1,:),fixtest(2,:));
axis([0 800 0 600]);
figure;scatter(fixInPlot(1,:),fixInPlot(2,:));
axis([0 800 0 600]);
figure;scatter(fixOutPlot(1,:),fixOutPlot(2,:));
axis([0 800 0 600]);
figure;scatter(regions(1,:),regions(2,:));
axis([0 800 0 600]);
% 
% error=[];
% for k=1:length(fixInPlot);
%     if fixInPlot(1,k)>200 && fixInPlot(1,k)<300
%         error=[error,k];
%     else
%     end
% end
errors=fixInPlot(:,error);
error_region=[];

fixInROI=zeros(2,length(fixtest));
for i=1:length(fixtest);
    for k=1:length(regions);
        if regions(1,k)==fixtest(1,i) && regions(2,k)==fixtest(2,i)
            fixInROI(:,i)=regions(:,k);
        else
        end
    end
end
fixInROI=fixInROI(fixInROI~=0);
fixInROI=reshape(fixInROI,2,(length(fixInROI)/2));

% save('SSCM_ROI_Fix_EyeDat_fixcell.mat','fixations','transitions');
%%

for subjloop=1:size(fixations,2);
    
end

%%


figure;plot(vel);hold on;
line([tempbeg(1) tempbeg(1) tempend(1) tempend(1)], [0 900 0 900],'LineStyle' , '--' , 'Color' , 'k','LineWidth' , 1);