%% Get Fixations for Each Trial

% Position of previous fixation > Position of current fixation > Position
% of next fixation 
% Time at start of fixation (ms), Time at end of fixation, Duration of fixation
% Length of inbound saccade, Length of outbound saccade
              % 1             2             3      4         5                6      
 fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
          'Inbound Time','Outbound Time','Duration',...
          'Inbound Saccade Length','Outbound Saccade Length',...
          'Fixation Number','Total Fixations','Transition In'};
          %     7               8             9
          %     10                         11
          %     12                 13              14

vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';


% Load list of stimulus conditions (1010:1:1189)
load([vardir 'SSCM90_Set1_cnd.mat']);
% Load X-Y Pixel Coordinates of All ROIs in Sets Specified Above
load([inddir 'SSCM_Set1.2.3.4.5.6.7_ROI_All_xy.mat']);
load([inddir 'SSCM_ROI_All_area.mat']);
% Load Eye Position During Stimulus Presentation for Sets Specified Above
% Recalibration with cp2tform affine transformation for files with 9-point
% cross and polynomial transformation for files with 63-point grid 
load([vardir 'SSCM_eyedat_12-Jul-2013_17h11m.mat']);
% Load List of Scene Numbers for ROIs
load([inddir 'SSCM_ROI_Scene.mat']);

% construct filter for eye velocity measure
fltord = 60;
lowpasfrq = 25;
nyqfrq = 1000 ./ 2;
flt = fir2(fltord,[0,lowpasfrq./nyqfrq,lowpasfrq./nyqfrq,1],[1,1,0,0]);

% Omit first fixation in center: fixStart=2, keep: fixStart=1
fixStart=2;


datFix=cell(size(eyedat));
for setloop=2:size(eyedat,1);
    
    datFix{setloop,1}=cell(1,size(eyedat{setloop,1},2));
    for subjloop=1:size(eyedat{setloop,1},2)
        
        datFix{setloop,1}=cell(1,size(eyedat{setloop,1},2));
        for fileloop=1:size(eyedat{setloop,1},2);
            
            datFix{setloop,1}{1,fileloop}=cell(size(eyedat{setloop,1}{1,fileloop},1),1);
            for trial=1:size(eyedat{setloop,1}{1,fileloop},1);
                scene=ROI_Scene{setloop,1}(trial,1);
                % Novel Presentation (+1010, +1011 for 2nd Presentation)
                cond=(((scene-1)*2)+1010);
                ind=find(cnd==cond);
                
                % differentiate and multiply with sampling rate to get velocity as deg/sec
                resampx=resample(eyedat{setloop,1}{1,fileloop}{trial}(1,:),5,1);
                resampy=resample(eyedat{setloop,1}{1,fileloop}{trial}(2,:),5,1);
                % calculate velocity
                x_v= diff(filtfilt(flt,1, resampx)) .* 1000;
                y_v= diff(filtfilt(flt,1, resampy)) .* 1000;
                % combine x- and y-velocity to get overall eye velocity
                vel = abs(complex(x_v,y_v));
                % lim = threshold for detecting saccade
                % Threshold of 19 captures short-distance saccades that often
                % fall on faces and connect distinct clusters of eye-movements
                lim = 30;
                % 30 is used by Klin group
                % Find beginning and end of saccades
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
                
                % GET:
                % Position of previous fixation > Position of current fixation > Position
                % of next fixation
                % Time at start of fixation (ms), Time at end of fixation, Duration of fixation
                % Length of inbound saccade, Length of outbound saccade
                for fixlop=fixStart:length(tempbeg)
                    % X coord
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3)=round((mean(resampx(tempbeg(fixlop):tempend(fixlop)))*24)+400); % fixation position, horizontal
                    % Y coord
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,4)=round((mean(resampy(tempbeg(fixlop):tempend(fixlop)))*24)+300); % fixation position, vertical
                    % Inbound time
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,7)=tempbeg(1,fixlop);
                    % Outbound time
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,8)=tempend(1,fixlop);
                    % Fixation duration
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,9)=datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,8)-datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,7);
                    % Fixation Number
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,12)=fixlop;
                    % Total # of Fixations
                    datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,13)=length(tempbeg);
                    % Check for transition out/in
                    % If first fixation, inbound fixation was from center cross
                    if fixlop==1
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1)=400;
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,2)=300;
                        % Calculate distance between center cross and 1st fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,10)=sqrt(sum(diff([datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1:2);datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3:4)],1,1).^2));
                    elseif fixlop~=1 && fixlop~=length(tempbeg)
                        % If 2nd or greater fixation, inbound position is
                        % position of previous fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1)=round((mean(resampx(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+400); % previous fixation position, horizontal
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,2)=round((mean(resampy(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+300); % previous fixation position, horizontal
                        % Outbound fixation for previous fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,5)=datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3);
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,6)=datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,4);
                        % Calculate inbound saccade length
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,10)=sqrt(sum(diff([datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1:2);datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3:4)],1,1).^2));
                        % Calculate outbound saccade length
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,11)=sqrt(sum(diff([datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,3:4);datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3:4)],1,1).^2));
                    elseif fixlop==length(tempbeg)
                        % If 2nd or greater fixation, inbound position is
                        % position of previous fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1)=round((mean(resampx(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+400); % previous fixation position, horizontal
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,2)=round((mean(resampy(tempbeg(fixlop-1):tempend(fixlop-1)))*24)+300); % previous fixation position, horizontal
                        % Outbound fixation for previous fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,5)=datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3);
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,6)=datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,4);
                        % Calculate inbound saccade length
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,10)=sqrt(sum(diff([datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,1:2);datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3:4)],1,1).^2));
                        % Calculate outbound saccade length
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,11)=sqrt(sum(diff([datFix{setloop,1}{1,fileloop}{trial,1}(fixlop-1,3:4);datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,3:4)],1,1).^2));
                        % Last fixation has no outbound fixation
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,5)=NaN;
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,6)=NaN;
                        % Last fixation has no outbound saccade length
                        datFix{setloop,1}{1,fileloop}{trial,1}(fixlop,11)=NaN;
                    end
                end
            end
        end
    end
end

dateSaved=date;
time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];
save([inddir 'SSCM_fix_',dateSaved,'_',time,'.mat'],'indOT','Sets','datFix');

%% Filter Fixations Through Set of ROIs: Whole-Item ROIs
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
% Load Fixation data
load([inddir 'SSCM_fix_17-Jun-2013_14h11m']);
% Load condition list 
load([vardir 'SSCM90_Set1_cnd.mat']);
% Load WholeItem ROI xy coordinates
load([inddir 'SSCM_Set1.2.3.4.5.6.7_ROI_All_xy.mat']);
% Load List of ROIs to include
load([inddir 'SSCM_ROI_Include.mat']);
% Load List of Scene #'s of ROIs
load([inddir 'SSCM_ROI_Scene.mat']);
load([inddir 'SSCM_indOT.mat']);

tic
fix_Whole=cell(size(datFix));
for setloop=2:size(datFix,1);
    
    fix_Whole{setloop,1}=cell(1,size(datFix{setloop,1},2));
    for fileloop=1:size(datFix{setloop,1},2);
        
        fix_Whole{setloop,1}{1,fileloop}=cell(size(ROI_All_xy{setloop,1},1),2);        
        % Novel Presentation (+1010, +1011 for 2nd Presentation)
        for p=1:2
            for roi=1:size(ROI_All_xy{setloop,1},1)
                scene=ROI_Scene{setloop,1}(roi,1);
                cond=(((scene-1)*2)+(1009+p));
                ind=find(cnd==cond);
                
                if ~isempty(ind)
                    regions=ROI_All_xy{setloop,1}{roi,1};
                    
                    r=1;
                    for fixlop=1:size(datFix{setloop,1}{1,fileloop}{ind,1},1)
                        if ismember(datFix{setloop,1}{1,fileloop}{ind,1}(fixlop,3:4),regions,'rows')
                            fix_Whole{setloop,1}{1,fileloop}{roi,p}(r,:)=datFix{setloop,1}{1,fileloop}{ind,1}(fixlop,:);
%                             disp(fixlop)
                            r=r+1;
                        else
                        end
                    end
                end
            end
        end
    end
end 
toc
load chirp 
sound(y,Fs)
clear y
clear Fs
save([inddir 'SSCM_fix_Whole.mat'],'fix_Whole','inddir','vardir','indOT'); 

%% Calculate transitions into Whole ROI
             % 1             2             3      4         5                6      
 fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
          'Inbound Time','Outbound Time','Duration',...
          'Inbound Saccade Length','Outbound Saccade Length',...
          'Fixation Number','Total Fixations','Transition In'};
          %     7               8             9
          %     10                         11
          %     12                 13              14

          numROI=1140;
tic
for setloop=2:size(fix_Whole,1);
    
    for fileloop=1:size(fix_Whole{setloop,1},2);
        
        for p=1:2
            for roi=1:numROI
                if ~isempty(fix_Whole{setloop,1}{1,fileloop}{roi,p})
                    fix_Whole{setloop,1}{1,fileloop}{roi,p}(:,14)=0;
                    
                    for fixlop=1:size(fix_Whole{setloop,1}{1,fileloop}{roi,p},1)
                        if fixlop>1
                            if diff(fix_Whole{setloop,1}{1,fileloop}{roi,p}(fixlop-1:fixlop,12))~=1
                                fix_Whole{setloop,1}{1,fileloop}{roi,p}(fixlop,14)=1;
                            else
                            end
                        else
                            fix_Whole{setloop,1}{1,fileloop}{roi,p}(fixlop,14)=1;
                        end
                    end
                end
            end
        end
    end
end 
toc

save([inddir 'SSCM_fix_Whole.mat'],'fix_Whole','fixLabl','inddir','vardir','indOT'); 

%% Create 'hit matrices' from fixation time/duration data: Whole-Item

%              1             2             3      4         5                6      
%  fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
%           'Inbound Time','Outbound Time','Duration',...
%           'Inbound Saccade Length','Outbound Saccade Length',...
%           'Fixation Number','Total Fixations'};
%               7               8             9
%               10                         11
%               12                 13
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
fixDir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\fixHit_Whole';

% load([inddir 'SSCM_fix_Whole.mat']);
% 
% for k=1:7
%     saveDir=[xydir num2str(k) '\'];   
%     ROI_All_areaSet=ROI_All_area{k,1};
%     save([saveDir 'SSCM_All_area.mat'],'ROI_All_areaSet');
% end

presTime=[10000,6000];

tic
fixHit_Whole=cell(size(fix_Whole));
for setloop=2:size(fix_Whole,1);
    saveDir=[fixDir num2str(setloop) '\']; 
    tic
    fixHit_Whole{setloop,1}=cell(1,size(fix_Whole{setloop,1},2));
    for fileloop=1:size(fix_Whole{setloop,1},2);
        
        fixHit_Whole{setloop,1}{1,fileloop}=cell(1140,2);
        % Novel Presentation (+1010, +1011 for 2nd Presentation)
        for p=1:2
            for roi=1:1140;
                fixHit_Whole{setloop,1}{1,fileloop}{roi,p}=zeros(1,presTime(p));
                for fixlop=1:size(fix_Whole{setloop,1}{1,fileloop}{roi,p},1)
                    % Start of fixation (ms from onset)
                    tStart=fix_Whole{setloop,1}{1,fileloop}{roi,p}(fixlop,7);
                    % End of fixation (ms from onset)
                    tEnd=fix_Whole{setloop,1}{1,fileloop}{roi,p}(fixlop,8);
                    if tEnd>presTime(p)
                        tEnd=presTime(p);
                    end
                    fixHit_Whole{setloop,1}{1,fileloop}{roi,p}(1,tStart:(tEnd-1))=1;
                end
            end
        end
    end
    if ~exist(saveDir,'dir');
        mkdir(saveDir);
    end
    save([saveDir 'SSCM_fixHit_Whole',num2str(setloop),'.mat'],'fixHit_Whole');
    clear fixHit_Whole
%     clear SSCM_All_area
    toc
end
toc
load chirp 
sound(y,Fs)
clear y
clear Fs
clear fix_Whole ROI_All_area ROI_All_xy 
% save([inddir 'SSCM_fixHit_Whole.mat'],'fixHit_Whole');



% Filter Fixations Through Set of ROIs: Face ROIs
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
% Load Fixation data
load([inddir 'SSCM_fix_17-Jun-2013_14h11m']);
% Load condition list 
load([vardir 'SSCM90_Set1_cnd.mat']);
% Load WholeItem ROI xy coordinates
load([inddir 'SSCM_ROI_Face_All_xy.mat']);
% Load List of ROIs to include
load([inddir 'SSCM_ROI_Include.mat']);
% Load List of Scene #'s of ROIs
load([inddir 'SSCM_ROI_Scene.mat']);
load([inddir 'SSCM_indOT.mat']);

tic
fix_Face=cell(size(datFix));
for setloop=2:size(datFix,1);
    
    fix_Face{setloop,1}=cell(1,size(datFix{setloop,1},2));
    for fileloop=1:size(datFix{setloop,1},2);
        
        fix_Face{setloop,1}{1,fileloop}=cell(size(ROI_Face_All_xy{setloop,1},1),2);        
        for p=1:2
            for roi=1:size(ROI_Face_All_xy{setloop,1},1)
                scene=ROI_Scene{setloop,1}(roi,1);
                cond=(((scene-1)*2)+(1009+p));
                ind=find(cnd==cond);
                
                if ~isempty(ind)
                    regions=ROI_Face_All_xy{setloop,1}{roi,1};
                    
                    r=1;
                    for fixlop=1:size(datFix{setloop,1}{1,fileloop}{ind,1},1)
                        if ismember(datFix{setloop,1}{1,fileloop}{ind,1}(fixlop,3:4),regions,'rows')
                            fix_Face{setloop,1}{1,fileloop}{roi,p}(r,:)=datFix{setloop,1}{1,fileloop}{ind,1}(fixlop,:);
%                             disp(fixlop)
                            r=r+1;
                        else
                        end
                    end
                end
            end
        end
    end
end 
toc
load chirp 
sound(y,Fs)
clear y
clear Fs
% save([inddir 'SSCM_fix_Face.mat'],'fix_Face','inddir');

% Calculate transitions into Face ROI
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
              % 1             2             3      4         5                6      
 fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
          'Inbound Time','Outbound Time','Duration',...
          'Inbound Saccade Length','Outbound Saccade Length',...
          'Fixation Number','Total Fixations','Transition In'};
          %     7               8             9
          %     10                         11
          %     12                 13              14
% 
% load([inddir 'SSCM_fix_Face.mat'])
load([inddir 'SSCM_indOT.mat']);

numROI=1140;
tic
for setloop=2:size(fix_Face,1);
    
    for fileloop=1:size(fix_Face{setloop,1},2);
        
        for p=1:2
            for roi=1:numROI
                if ~isempty(fix_Face{setloop,1}{1,fileloop}{roi,p})
                    
                    fix_Face{setloop,1}{1,fileloop}{roi,p}(:,14)=0;
                    for fixlop=1:size(fix_Face{setloop,1}{1,fileloop}{roi,p},1)
                        if fixlop>1
                            if diff(fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop-1:fixlop,12))~=1
                                fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,14)=1;
                            else
                            end
                        else
                            fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,14)=1;
                        end
                    end
                end
            end
        end
    end
end 
toc

save([inddir 'SSCM_fix_Face.mat'],'fix_Face','fixLabl','inddir','vardir','indOT'); 

% Create 'hit matrices' from fixation time\duration data: Face Only

%              1             2             3      4         5                6      
%  fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
%           'Inbound Time','Outbound Time','Duration',...
%           'Inbound Saccade Length','Outbound Saccade Length',...
%           'Fixation Number','Total Fixations'};
%               7               8             9
%               10                         11
%               12                 13
vardir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\';
inddir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\';
fixDir='S:\Drew\Backup_7.1.13\Buffalo Rotation\Scene Manipulation\Variables\ROI_Ind\Files\fixHit_Face';

% load([inddir 'SSCM_fix_Face.mat']);
% 
% for k=1:7
%     saveDir=[xydir num2str(k) '\'];   
%     ROI_All_areaSet=ROI_All_area{k,1};
%     save([saveDir 'SSCM_All_area.mat'],'ROI_All_areaSet');
% end

presTime=[10000,6000];

tic
fixHit_Whole=cell(size(fix_Face));
for setloop=2:size(fix_Face,1);
    saveDir=[fixDir num2str(setloop) '\']; 
    tic
    fixHit_Face{setloop,1}=cell(1,size(fix_Face{setloop,1},2));
    for fileloop=1:size(fix_Face{setloop,1},2);
        
        fixHit_Face{setloop,1}{1,fileloop}=cell(1140,2);
        % Novel Presentation (+1010, +1011 for 2nd Presentation)
        for p=1:2
            for roi=1:1140;
                fixHit_Face{setloop,1}{1,fileloop}{roi,p}=zeros(1,presTime(p));
                for fixlop=1:size(fix_Face{setloop,1}{1,fileloop}{roi,p},1)
                    % Start of fixation (ms from onset)
                    tStart=fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,7);
                    % End of fixation (ms from onset)
                    tEnd=fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,8);
                    if tEnd>presTime(p)
                        tEnd=presTime(p);
                    end
                    fixHit_Face{setloop,1}{1,fileloop}{roi,p}(1,tStart:(tEnd-1))=1;
                end
            end
        end
    end
    if ~exist(saveDir,'dir');
        mkdir(saveDir);
    end
    save([saveDir 'SSCM_fixHit_Face',num2str(setloop),'.mat'],'fixHit_Face');
    clear fixHit_Face
%     clear SSCM_All_area
    toc
end
toc
load chirp 
sound(y,Fs)
clear y
clear Fs

% save([inddir 'SSCM_fixHit_Whole.mat'],'fixHit_Whole');





%% Create 'hit matrices' from fixation time\duration data: Faces

%              1             2             3      4         5                6      
%  fixLabl={'Inbound xPos','Inbound yPos','xPos','yPos','Outbound xPos','Outbound yPos',...                   
%           'Inbound Time','Outbound Time','Duration',...
%           'Inbound Saccade Length','Outbound Saccade Length',...
%           'Fixation Number','Total Fixations'};
%               7               8             9
%               10                         11
%               12                 13
presTime=[10000,6000];

fixHit_Face=cell(size(fix_Face));
for setloop=2:size(fix_Face,1);
    
    fixHit_Face{setloop,1}=cell(1,size(fix_Face{setloop,1},2));
    for fileloop=1:size(fix_Face{setloop,1},2);
        
        fixHit_Face{setloop,1}{1,fileloop}=cell(size(ROI_Face_All_xy{setloop,1},1),2);
        % Novel Presentation (+1010, +1011 for 2nd Presentation)
        for p=1:2
            for roi=1:size(ROI_Face_All_xy{setloop,1},1) 
                fixHit_Face{setloop,1}{1,fileloop}{roi,p}=zeros(1,presTime(p));
                for fixlop=1:size(fix_Face{setloop,1}{1,fileloop}{roi,p},1)
                    % Start of fixation (ms from onset)
                    tStart=fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,7);
                    % End of fixation (ms from onset)
                    tEnd=fix_Face{setloop,1}{1,fileloop}{roi,p}(fixlop,8);
                    fixHit_Face{setloop,1}{1,fileloop}{roi,p}(1,tStart:(tEnd-1))=1;
                end
            end
        end
    end
end

save([inddir 'SSCM_fixHit_Face.mat'],'fixHit_Face','indOT','inddir','vardir');

