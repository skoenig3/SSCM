%% Import Face ROIs and Get Scene Numbers
% Get List of Folders with Face ROIs for each Set of Scenes
% ROI_FaceDir=getfiles('/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/Face ROI Scene/');
% Get List of Folders with Whole-Head ROIs for each Set of Scenes
faceDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/WholeHeadROIs/';
ROI_FaceDir=getfiles(faceDir);

% Return only folders with Face ROIs
ROI_FaceDir=ROI_FaceDir(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_FaceDir, 'FaceROI')));

% Get Face ROIs from each Set folder 
ROI_Face=cell(length(ROI_FaceDir),1);
for k=1:length(ROI_FaceDir)
    x=getfiles([faceDir,ROI_FaceDir{k,1},'/']);
    ROI_Face{k,1}=x(cellfun(@(IDX) ~isempty(IDX), strfind(x, 'ROI')));
    clear x
end

% Find Scene Number for Each Face ROI in Each Set
for k=1:length(ROI_Face);
    for j=1:length(ROI_Face{k,1});
        ROI_FaceSNum{k,1}(j,1)=str2double(ROI_Face{k,1}{j,1}(13:14));
    end
end
%% Get X-Y Coordinates of Pixels in Each Face ROIs
ROI_Face_xy=cell(length(ROI_Face),1);
ROI_Face_area=cell(length(ROI_Face),1);
% For each Set of Scenes...
tic
for setloop=1:length(ROI_Face);
    tic
    % Ignore Sets without Face ROIS (Set1 as of 12.12.04)
    if ~isempty(ROI_Face{setloop,1})
        % Find X-Y Coordinates of Pixels in Each Face ROI for Set setloop
        for k=1:length(ROI_Face{setloop,1});
            % Concatenate String of Current Set Number with FaceROI tag
            setStr=[ROI_FaceDir{setloop,1}(2:3),'_FaceROI'];
            % Read black and white .bmp of k Face ROI for Set setloop
            I = imread([faceDir,'S',setStr,'/',ROI_Face{setloop,1}{k,1}]);
            % Convert black and white .bmp of Face ROI to grayscale image
            Igray=rgb2gray(I);
            % Clear image I to %%save memory
            clear I
            % Get X-Y coordinates (row & column) of pixels with white (member of ROI)
            [Y,X]=find(Igray > 0);
            % Cleary Igray to %%save memory
            clear Igray
            % Store X-Y Coordinates After Subtracting 600 pixels from Y
            % coordinates because ROI is flipped 
            ROI_Face_xy{setloop,1}{k,1}=[X 600-Y];
            % Store Total Area Occupied by Face ROI (# of pixels)
            ROI_Face_area{setloop,1}(k,1) = length(ROI_Face_xy{setloop,1}{k,1});
            clear X
            clear Y
        end
    else
    end
    toc
end
toc
%% Separately calculate Face ROI Area
ROI_Face_area=cell(length(ROI_Face_All_xy),1);
for setloop=1:length(ROI_Face_All_xy);
    % Ignore Sets without Face ROIS (Set1 as of 12.12.04)
    if ~isempty(ROI_Face_All_xy{setloop,1})
        % Find X-Y Coordinates of Pixels in Each Face ROI for Set setloop
        for k=1:length(ROI_Face_All_xy{setloop,1});
            % Store Total Area Occupied by Face ROI (# of pixels)
            ROI_Face_area{setloop,1}(k,1) = length(ROI_Face_All_xy{setloop,1}{k,1});
        end
    else
    end
end

%%
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';

date%%saved=date;
time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];
%%save([inddir 'SSCM_ROI_Face_area_',date%%saved,'_',time,'.mat'],'ROI_Face_area','Sets','date%%saved','time');
%%
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';

ROI_All_area=cell(length(ROI_Face_All_xy),1);
for setloop=1:length(ROI_All_xy);
    % Ignore Sets without Face ROIS (Set1 as of 12.12.04)
    if ~isempty(ROI_All_xy{setloop,1})
        % Find X-Y Coordinates of Pixels in Each Face ROI for Set setloop
        for k=1:length(ROI_All_xy{setloop,1});
            % Store Total Area Occupied by Face ROI (# of pixels)
            ROI_All_area{setloop,1}(k,1) = length(ROI_All_xy{setloop,1}{k,1});
        end
    else
    end
end
%%save([inddir 'SSCM_ROI_All_area.mat'],'ROI_All_area','Sets','date%%saved','time');

%% Load ROIs and XY pixel coordinates
varDir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/';
roidir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/Set';
SetsNum=[1,2,3,4,5,6,7];

for setloop=1:length(SetsNum)
    if SetsNum(setloop)<10
        setString=['0' num2str(SetsNum(setloop))];
    else
        setString=num2str(SetsNum(setloop));
    end
    ROI_All{setloop,1}=getfiles([roidir,setString,'/'],'*.bmp');
end

% Fruit ROIs
% fruit{setloop,1}(f,1)=strfind(ROI_All{setloop,1}{1,fileloop}{roi,1},'fruit');
% fruit{setloop,1}(f,1)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{setloop,1}{1,fileloop}, 'fruit'));
% indTrl{Set,1}(:,1)=cellfun(@(IDX) isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACEWITH_'));

for setloop=2:size(ROI_All,1);
    fruitInd{setloop,1}=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{setloop,1}, 'fruit'));
    fruitStr{setloop,1}=ROI_All{setloop,1}(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{setloop,1}, 'fruit')));
    disp(['Set ',num2str(setloop),': ',num2str(sum(fruitInd{setloop,1})),' fruit'])
end
%%save([varDir,'SSCM_fruitInd_S1-7.mat'],'fruitInd');
%%
% Place Face ROIs into same position as source Monkey and empty cell
% otherwise so that the same logical indices can be used on the Face ROIs
% (ROI_Cat for item category and indTrl for ROI trial type)

tic
ROI_Face_All_xy=cell(length(ROI_Face),1);
ROI_Face_All_area=cell(length(ROI_Face),1);
% Go through every set of Face ROIs
% SKIP SET 1
for setloop=2:length(ROI_All);
    % Loop through every whole item ROI for a given set
    ROI_Face_All_xy{setloop,1}=cell(length(ROI_All{setloop,1}),1);
    ROI_Face_All_area{setloop,1}=NaN(length(ROI_All{setloop,1}),1);
    for allLoop=1:length(ROI_All{setloop,1});
        % Create string from Whole Item ROI that matches Face ROI
        faceStr=[ROI_All{setloop,1}{allLoop,1}(1:18),'ROI_Face_',ROI_All{setloop,1}{allLoop,1}(19:end)];
        % Replace 'REPLACEWITH_' Tag to Match Face ROI
        if ~isempty(strfind(faceStr,'REPLACEWITH_'));
            faceStr=strrep(faceStr,'REPLACEWITH_','');
        % Replace 'REPLACE_' Tag to Match Face ROI
        elseif ~isempty(strfind(faceStr,'REPLACE_'));
            faceStr=strrep(faceStr,'REPLACE_','');
        % Replace 'REPEAT_' Tag to Match Face ROI
        elseif ~isempty(strfind(faceStr,'REPEAT_'));
            faceStr=strrep(faceStr,'REPEAT_','');
        else
        end
        % Find Location of Face ROI that matches Whole Item ROI
        ind=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_Face{setloop,1},faceStr));
        % If the Whole Item ROI is a Face ROI, add x-y coordinates to
        % position of Whole Item ROI
        if sum(ind)==1;
            ROI_Face_All_xy{setloop,1}(allLoop,1)=ROI_Face_xy{setloop,1}(ind,1);
            ROI_Face_All_area{setloop,1}(allLoop,1)=ROI_Face_area{setloop,1}(ind,1);
        else
        end
    end
end
toc

%% %%save xy coordinates for Face ROIs
date%%saved=date;
time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];

inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';


%%save([inddir 'SSCM_ROI_Face_All_xy_',date%%saved,'_',time,'.mat'],'ROI_Face_All_xy','Sets','SetsNum');
%%save([inddir 'SSCM_ROI_Face_All_area_',date%%saved,'_',time,'.mat'],'ROI_Face_All_area','Sets','SetsNum');

%%
mNum=cell(length(ROI_Face_All_xy),1);
for k=1:length(ROI_Face_All_xy);
    if ~isempty(ROI_Face_All_xy{k,1});
        mNum{k,1}=ROI_Face_All_xy{k,1}(ROI_Cat{k,1}{2,1}(:,1) & indTrl{k,1}(:,1));
        indNum(k,1)=sum((ROI_Cat{k,1}{2,1}(:,1) & indTrl{k,1}(:,1)),1);
        indNum(k,2)=sum(indTrl{k,1}(:,1),1);
        indNum(k,3)=sum(ROI_Cat{k,1}{2,1}(:,1),1);
%         mNum{k,1}=ROI_Face_All_xy{k,1}(indTrl{k,1}(:,1));
    else
    end
end

for setloop=2:length(ROI_Face_area);
    faceMean(setloop,1)=nanmean(ROI_Face_area{setloop,1});
end


%% Import ROI images, Sort them and Get xy pixel coordinates
%What set?
Set=3;
vardir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/';
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';
% Load previous Sets of ROIs, ROI coordinates and trial types
matlist=getfiles(inddir,'*.mat');
for k=1:length(matlist)
    load([inddir matlist{k}]);
end
SetsNum(1,Set)=Set;


%%
SceneFolders=getfiles(['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/Scenes/Set',num2str(Set)]);
RepT=SceneFolders(cellfun(@(IDX) ~isempty(IDX), strfind(SceneFolders, 'Repeated')));
MR=SceneFolders(cellfun(@(IDX) ~isempty(IDX), strfind(SceneFolders, 'MonkeyReplaced')));
OR=SceneFolders(cellfun(@(IDX) ~isempty(IDX), strfind(SceneFolders, 'ObjectReplaced')));
for k=1:length(MR)
    trlType{Set,1}(k,1)=str2double(RepT{k,1}(6:7));
    trlType{Set,1}(k,2)=str2double(MR{k,1}(6:7));
    trlType{Set,1}(k,3)=str2double(OR{k,1}(6:7));
end
%%
% Identify Scene Number for each ROI
for Set=1:length(SetsNum);    
    for k=1:length(ROI_All{Set,1});
        ROI_Scene{Set,1}(k,1) = str2double(ROI_All{Set,1}{k}(13:14));
    end
end

for Set=1:length(SetsNum);
    % Sort Monkeys and Objects
    ROI_Item_Labl={'Monkeys','Objects'};
    ROI_Item = zeros(length(ROI_All{Set,1}),length(ROI_Item_Labl));
    for k=1:length(ROI_All{Set,1});
        %Match Un-manipulated ROIs    // Match for Replace ROIs     // Match ReplaceWith          // Match Repeat
        if ROI_All{Set,1}{k}(19)=='M' || ROI_All{Set,1}{k}(27)=='M' || ROI_All{Set,1}{k}(31)=='M' || ROI_All{Set,1}{k}(26)=='M'
            ROI_Item(k,1) = 1;
        elseif ROI_All{Set,1}{k}(19)=='O' || ROI_All{Set,1}{k}(27)=='O' || ROI_All{Set,1}{k}(31)=='O' || ROI_All{Set,1}{k}(26)=='O'
            ROI_Item(k,2) = 1;
        end
    end
    ROI_Cat{Set,1}{2,1}=logical(ROI_Item);
    clear ROI_Item
end

ROI_Cat_Labl=ROI_Cat_Labl{1,1};
ROI_Cat_Labl{6,1}='Gaze';
ROI_Cat_Labl{6,2}={'Direct','Averted'};

%%


% Index of All ROIs shown in 1st Presentation
indTrl{Set,1}(:,1)=cellfun(@(IDX) isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACEWITH_'));
% nov=ROI_Scene(ind.trl.Nov);
% Index of All ROIs shown in 1st Presentation
indTrl{Set,1}(:,2)=cellfun(@(IDX) isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACE_'));
% rep=ROI_Scene(ind.trl.Rep);
for k=1:length(ROI_All{Set,1})
    % Index of All Repeated w/o Manipulation Scene ROIs
    if ismember(str2double(ROI_All{Set,1}{k}(13:14)),trlType{Set,1}(:,1));
        indTrl{Set,1}(k,3)=true;
        indTrl{Set,1}(k,4)=false;
        indTrl{Set,1}(k,5)=false;
    % Index of All Monkey Replaced Scene ROIs
    elseif ismember(str2double(ROI_All{Set,1}{k}(13:14)),trlType{Set,1}(:,2));
        indTrl{Set,1}(k,3)=false;
        indTrl{Set,1}(k,4)=true;
        indTrl{Set,1}(k,5)=false;
    % Index of All Object Replaced Scene ROIs
    elseif ismember(str2double(ROI_All{Set,1}{k}(13:14)),trlType{Set,1}(:,3));
        indTrl{Set,1}(k,3)=false;
        indTrl{Set,1}(k,4)=false;
        indTrl{Set,1}(k,5)=true;
    else
        disp('No match for trial type')
    end
end
% Index of Replaced ROI
% Replace M: CR1
indTrl{Set,1}(:,6)=(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACE_')) & ROI_Item(:,1));
% Replace O: CR1
indTrl{Set,1}(:,7)=(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACE_')) & ROI_Item(:,2));
% Replace M: CR2
indTrl{Set,1}(:,8)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACEWITH_')) & ROI_Item(:,1);
% Replace O: CR2
indTrl{Set,1}(:,9)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPLACEWITH_')) & ROI_Item(:,2);



% Sort by ROI Area (1=10000 pixels, 2=5000, 3=2000)
ROI_Size_Labl={'10000 pixels','5000 pixels','2000 pixels'};
ROI_Size(:,1)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '__1.'));
ROI_Size(:,2)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '__2.'));
ROI_Size(:,3)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '__3.'));




% Find Monkeys with 2, 1 or 0 Eyes Visible
ROI_M_Eye_Labl = {'No Eyes Visible','1 Eye Visible','2 Eye Visible'};
ROI_M_Eye(:,1)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_0_'));
ROI_M_Eye(:,2)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_1_'));
ROI_M_Eye(:,3)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_2_'));





% Separate ROIs by Age: Logical arrays
ROI_M_Age_Labl = {'Infant','Juvenile','Adult'};
ROI_M_Age(:,1)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_I_'));
ROI_M_Age(:,2)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_J_'));
ROI_M_Age(:,3)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_A_'));




%Separate by Sex (as surmised by Drew Solyst...but not confirmed by anyone)
ROI_M_Sex_Labl = {'Male','Female','Unknown'};
ROI_M_Sex(:,1)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_M.'));
ROI_M_Sex(:,2)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_F.'));
ROI_M_Sex(:,3)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, '_U.'));



% Find Monkeys with Direct or Averted Gaze
for Set=1:length(SetsNum);   
    ROI_Cat{Set,1}{6,1}=[];
    for k=1:length(ROI_All{Set,1});
           %ROI is 10000 pixels   and is a Monkey
        if (ROI_Cat{Set,1}{1,1}(k,1) && ROI_Cat{Set,1}{2,1}(k,1))
            if strcmp(ROI_All{Set,1}{k,1}(37:38),'Di');
                ROI_Cat{Set,1}{6,1}(k,1)=true;
                ROI_Cat{Set,1}{6,1}(k,2)=false;
            elseif ~strcmp(ROI_All{Set,1}{k,1}(37),'D')
                ROI_Cat{Set,1}{6,1}(k,1)=false;
                ROI_Cat{Set,1}{6,1}(k,2)=true;
            else
            end
        else
            ROI_Cat{Set,1}{6,1}(k,1)=false;
            ROI_Cat{Set,1}{6,1}(k,2)=false;
        end
    end
end

% Find Repeat w/o Manipulation ROIs
indTrl{Set,1}(:,10)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPEAT_M'));
indTrl{Set,1}(:,11)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPEAT_O'));



ROI_Cat{Set,1}=[{ROI_Size};{ROI_Item};{ROI_M_Eye};{ROI_M_Age};{ROI_M_Sex};{ROI_Gaze}];


%% Get X-Y Coordinates of Pixels in All ROIs
roidir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/Set';
load([vardir 'ROI_All_Set1_xy_labl.mat']);
ROI_All_xy={ROI_All_xy};
ROI_All_area={ROI_All_area};
for setloop=1:length(ROI_All);
    for k=1:length(ROI_All{setloop,1});
        if SetsNum(setloop)<10
            setString=['0' num2str(SetsNum(setloop))];
        else
            setString=num2str(SetsNum(setloop));
        end
        I = imread([roidir,setString,'/',ROI_All{SetsNum(setloop),1}{k,1}]);
        Igray=rgb2gray(I);
        clear I
        [Y,X]=find(Igray > 0);
        clear Igray
        ROI_All_xy{setloop,1}{k,1}=[X 600-Y];
        ROI_All_area{setloop,1}(k,1) = length(ROI_All_xy{setloop,1}{k,1});
        clear X
        clear Y
    end
end

%%

Setstr='Set';
for k=1:length(SetsNum);
    if k<length(SetsNum);
        Setstr=[Setstr,num2str(SetsNum(k)),'.'];
    elseif k==length(SetsNum);
        Setstr=[Setstr,num2str(SetsNum(k))];
    end
end
date%%saved=date;

time=clock;
time=[num2str(time(4)),'h',num2str(time(5)),'m'];
%%
copydir=['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind_Copy/',date%%saved,'_',time];
if ~exist(copydir,'dir');
    mkdir(copydir);
end
movefile(inddir,copydir);

if ~exist(inddir,'dir');
    mkdir(inddir);
end

%%save([inddir 'SSCM_',Setstr,'_ROI_All_xy.mat'],'ROI_All_area','ROI_All_xy','SetsNum','date%%saved');
%%save([inddir 'SSCM_',Setstr,'_indTrl.mat'],'indTrl','indTrl_Labl','trlType','trlType_Labl','SetsNum','date%%saved');
%%save([inddir 'SSCM_',Setstr,'_ROI_Cat.mat'],'ROI_Cat','ROI_Cat_Labl','trlType','trlType_Labl','SetsNum','date%%saved','time');
%%save([inddir 'SSCM_',Setstr,'_indTrl_ROI_Cat.mat'],'ROI_Scene','ROI_Cat','ROI_Cat_Labl','indTrl','indTrl_Labl','trlType','trlType_Labl','SetsNum','date%%saved','time');
%%save([inddir 'SSCM_',Setstr,'_ROI_All.mat'],'ROI_All','SetsNum','date%%saved','time');

%% %%save Face ROIs
inddir='/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Variables/ROI_Ind/Files/';

%%save([inddir,'SSCM_',Setstr,'ROI_Face_xy.mat'],'ROI_Face_xy'

%% Get number of ROIs indexed for each category
catNum=cell(6,1);
for cat=1:length(ROI_Cat{setloop,1})
    for level=1:size(ROI_Cat{1,1}{cat,1},2);
        for setloop=1:length(ROI_Cat)
            catNum{cat,1}(setloop,level)=sum(ROI_Cat{setloop,1}{cat,1}(:,level));
        end
    end
end

%% Add Direct & Averted Gaze Indices to ROI_Cat

ROI_Cat_LablT=ROI_Cat_Labl{1,1};
ROI_Cat_LablT{6,1}='Gaze';
ROI_Cat_LablT{6,2}={'Direct','Averted'};
ROI_Cat{1,1}{6,1}=ROI_Gaze{1,1};
ROI_Cat{2,1}{6,1}=ROI_Gaze{2,1};

%% Plot Histogram of Face ROI Area
ROIarea=[];
for k=1:length(ROI_Face_area)
    ROIarea=[ROIarea;ROI_Face_area{k,1}];
end
figure;
hist(ROIarea);
nanmean(ROIarea)
%% Adjust to fit in more sets
ROI_Cat={ROI_Cat};
ROI_Cat_Labl={ROI_Cat_Labl};
SetsNum=[1,2];
indTrl{Set,1}(:,10)=ROI_M_RepT_CR;
indTrl{Set,1}(:,11)=ROI_O_RepT_CR;
indTrl={indTrl};
indTrl_Labl={'First Pres.','Second Pres.','Repeat','Replace M','Replace O','Replace M: CR1','Replace O: CR1'...
    'Replace M: CR2','Replace O: CR2','Repeat M CR','Repeat O CR'};
ROI_All{Set,1}=getfiles('/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/Set1','*.bmp');
% trlType{1,1}={trlType{1,1},{trlType{1,2},{trlType{1,3}};
type=cell(size(SetsNum));
type{1,1}={trlType{1,1},trlType{1,2},trlType{1,3}};
type{2,1}={trlType{2,1},trlType{2,2},trlType{2,3}};
trlType=type;
trlType_Labl={'Repeat','Monkey Replaced','Object Replaced'};
%%save([vardir 'S1.S2_trltype_11.09.12.mat'],'trlType','trlType_Labl');
%% Load ROIs and XY pixel coordinates
ROI_All{Set,1}=getfiles(['/Users/DrewSolyst/Documents/Buffalo Rotation/Scene Manipulation/Stimuli/ROIs/Set',num2str(Set)],'*.bmp');

%% Sort out Set 1 trlTyp
Monkrepeat=Monkrepeat(4,:);
Objrepeat=Objrepeat(4,:);
Monkrepl=Monkrepl(4,:);
Objrepl=Objrepl(4,:);

load([vardir 'CreateScene_S',num2str(Set),'_trltyp.mat']);
sRepeat=sort([Monkrepeat';Objrepeat']);
sRepeatStr=cell(size(sRepeat));
for k=1:length(sRepeat)
    if sRepeat(k)<10
        sRepeatStr{k,1}=['Scene0',num2str(sRepeat(k))];
    else
        sRepeatStr{k,1}=['Scene',num2str(sRepeat(k))];
    end
end
sRepeatStr=num2str(sRepeat);
sReplace=sort([Monkrepl';Objrepl']);
sMRepl=sort(Monkrepl)';
sORepl=sort(Objrepl)';
trlType={[sRepeat,sMRepl,sORepl]};
trlType1=trlType;
trlType2={trlType1{1,1}(:,1),trlType1{1,1}(:,2),trlType1{1,1}(:,3)};
trlTypex=cell(2,3);
trlTypex(2,:)=trlType2;

trlTypex(1,:)={RepT,MR,OR};
trlType=trlTypex;
%%save([vardir 'S1.S2_trltype_11.09.12.mat'],'trlType','trlType_Labl');

typex=cell(2,1);
typex{1,1}=[trlType{1,1}{1,1},trlType{1,1}{1,2},trlType{1,1}{1,3}];
typex{2,1}=[trlType{2,1}{1,1},trlType{2,1}{1,2},trlType{2,1}{1,3}];
trlType=typex;


%% Sets 1, 2 & 3 did not properly label Repeat w/o manipulation ROIs
% Sets 1 & 2 had no label, Set 3 only attached Repeat label to M's in MT
% and O's in OT scenes.
% 
% sRepeatStr=num2str(trlType{Set,1}(k,1));
% raandom=randi(2,[47,1]);
% % Select Repeated w/o manipulation Monkeys (2 eyes visible, 5000 pixel Monkeys)
% ROI_M_RepT_2=ROI_All{Set,1}(indTrl{Set,1}(:,3) & ROI_M_Eye(:,3) & ROI_Size(:,2));
% ROI_M_RepT_CRstr=cell(size(trlType{Set,1}(:,1)));
% for k=1:length(trlType{Set,1}(:,1));
%     dum=cell(2,1);
%     z=1;
%     for j=1:length(ROI_M_RepT_2)
%         if str2double(ROI_M_RepT_2{j,1}(13:14))==trlType{Set,1}(k,1);
%            dum{z,1}=ROI_M_RepT_2{j,1};
%            z=z+1;
%         else
%         end
%     end
%     if ~isempty(dum{2,1});
%         x=randi(2,1);
%         ROI_M_RepT_CRstr{k,1}=dum{x,1};
%     else
%         ROI_M_RepT_CRstr{k,1}=dum{1,1};
%     end
% end
% ROI_M_RepT_CR=ismember(ROI_All{Set,1},ROI_M_RepT_CRstr);
% indTrl{2,1}(:,10)=ROI_M_RepT_CR;
% 
% ROI_O_RepT=ROI_All{2,1}(indTrl{2,1}(:,3) & ROI_Item(:,2) & ROI_Size(:,2));
% ROI_O_RepTsel=ROI_O_RepT(1:2:length(ROI_O_RepT),1);
% ROI_O_RepT_CR=ismember(ROI_All{Set,1},ROI_O_RepTsel);
% indTrl{2,1}(:,11)=ROI_O_RepT_CR;
% 
% %Start random counter
% sRepeatStr=num2str(trlType{Set,1}(k,1));
% raandom=randi(2,[47,1]);
% % Select Repeated w/o manipulation Monkeys (2 eyes visible, 5000 pixel Monkeys)
% ROI_M_RepT_2=ROI_All{Set,1}(indTrl{Set,1}(:,3) & ROI_M_Eye(:,3) & ROI_Size(:,2));
% ROI_M_RepT_CRstr=cell(size(trlType{Set,1}(:,1)));
% for k=1:length(trlType{Set,1}(:,1));
%     dum=cell(2,1);
%     z=1;
%     for j=1:length(ROI_M_RepT_2)
%         x=strcmp(ROI_M_RepT_CRstr,ROI_M_RepT_2{j,1});
%         if sum(x)<1
%             ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_2{j,1};
%         elseif sum(x)>0
%             if raandom(j)==1
%                 ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_2{j,1};
%             else
%             end
%             
%         else
%         end
%     end
% end
% 
% ROI_M_RepT_CR=ismember(ROI_All{Set,1},ROI_M_RepT_CRstr);
% indTrl(:,10)=ROI_M_RepT_CR;

% Select Repeat ROIs
indTrl{Set,1}(:,10)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPEAT_M'));
indTrl{Set,1}(:,11)=cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{Set,1}, 'REPEAT_O'));
% Set 3 Added "REPEAT_" only to 15 objects and 15 monkeys    
ROI_M_RepT_2=ROI_All{Set,1}(indTrl{Set,1}(:,3) & ROI_M_Eye(:,3) & ROI_Size(:,2));
ROI_M_RepT_T=ROI_M_RepT_2(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_M_RepT_2, 'REPEAT_M')));
ROI_M_RepT_T2=ROI_M_RepT_2(cellfun(@(IDX) isempty(IDX), strfind(ROI_M_RepT_2, 'REPEAT_M')));
ROI_M_RepT_CRstr=cell(15,1);
for k=1:15;
    for j=1:length(ROI_M_RepT_T);
        if k>1
            x=strcmp(ROI_M_RepT_CRstr{k-1}(13:14),ROI_M_RepT_T{j,1}(13:14));
            if sum(x)<1
                ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_T{j,1};
            end
        else
            ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_2{j,1};
        end
    end
end
ROI_M_RepT_Sel=ROI_M_RepT_T2([3,4,6,9,10,14,15,18,20,21,22,27,30],1);
for k=1:length(ROI_M_RepT_T)
    RepeatNum(k,1)=str2double(ROI_M_RepT_T{k,1}(13:14));
end

ROI_M_RepT_Sel=cell(length(ROI_M_RepT_T2),1);
z=1;
for k=1:length(ROI_M_RepT_T2)
    if ~ismember(str2double(ROI_M_RepT_T2{k,1}(13:14)),RepeatNum)
        ROI_M_RepT_Sel{z,1}=ROI_M_RepT_T2{k,1};
        z=z+1;
    end
end


% Select Repeated w/o manipulation Monkeys (2 eyes visible, 5000 pixel Monkeys)
ROI_M_RepT_2=ROI_All{Set,1}(indTrl{Set,1}(:,3) & ROI_M_Eye(:,3) & ROI_Size(:,2));
ROI_M_RepT_CRstr=cell(size(trlType{Set,1}(:,1)));
for k=1:length(trlType{Set,1}(:,1));
    dum=cell(2,1);
    z=1;
    for j=1:length(ROI_M_RepT_2)
        x=strcmp(ROI_M_RepT_CRstr,ROI_M_RepT_2{j,1});
        if sum(x)<1
            ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_2{j,1};
        elseif sum(x)>0
            if raandom(j)==1
                ROI_M_RepT_CRstr{k,1}=ROI_M_RepT_2{j,1};
            else
            end
            
        else
        end
    end
end

RepeatM=cell(35,1);
for k=1:length(ROI_M_RepT_T2)
    if ~ismember(str2double(ROI_M_RepT_T2{k,1}(13:14)),RepeatNum);
        RepeatM{k,1}=ROI_M_RepT_T2{k,1}(:,1:end);
    else
    end
end
% Manually chosen Repeat M ROIs (2.M_2 eyes in Repeat scenes with lowest ID #)
x=[3,4,6,9,10,14,15,18,20,21,22,27,30,32,33];
RepeatM2=ROI_M_RepT_T2([3,4,6,9,10,14,15,18,20,21,22,27,30,32,33],1);
RepeatM=[RepeatM2;ROI_M_RepT_T];
ROI_M_RepT_CR=ismember(ROI_All{Set,1},RepeatM);
indTrl{3,1}(:,10)=ROI_M_RepT_CR;

% Choose Repeat O's
ROI_O_RepT=ROI_All{3,1}(cellfun(@(IDX) ~isempty(IDX), strfind(ROI_All{3,1}, 'REPEAT_O')));
for k=1:length(ROI_M_RepT_T)
    RepeatNumO(k,1)=str2double(ROI_O_RepT{k,1}(13:14));
end

RepeatO=cell(length(ROI_O_RepT2),1);
for k=1:length(ROI_O_RepT2)
    if ~ismember(str2double(ROI_O_RepT2{k,1}(13:14)),RepeatNumO);
        RepeatO{k,1}=ROI_O_RepT2{k,1}(:,1:end);
    else
    end
end
% Manually chosen Repeat O ROIs (2.O in Repeat scenes with lowest ID #)
y=[1,3,8,12,14,16,19,22,25,30,32,34,37,42,44];
RepeatO=[ROI_O_RepT2([1,3,8,12,14,16,19,22,25,30,32,34,37,42,44],1);ROI_O_RepT];
% ROI_O_RepT2=ROI_All{3,1}(indTrl{3,1}(:,3) & ROI_Item(:,2) & ROI_Size(:,2));
% ROI_O_RepTsel=ROI_O_RepT(1:2:length(ROI_O_RepT),1);
ROI_O_RepT_CR=ismember(ROI_All{Set,1},RepeatO);
indTrl{3,1}(:,11)=ROI_O_RepT_CR;
   
